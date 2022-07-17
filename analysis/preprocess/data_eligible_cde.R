################################################################################

# This script:
# - reads the processed data
# - applies eligibility criteria from boxes c and d of Figure 3 in protocol
# - saves processed data from eligible individuals

################################################################################

# import libraries ----
library(tidyverse)
library(lubridate)
library(glue)

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

# individuals who are eligible based on criteria in box a of Figure 3 on protocol
data_eligible_a <- readr::read_rds(
  here::here("output", "data", "data_eligible_a.rds")) 

# individuals who are eligible based on criteria in box b of Figure 3 on protocol
data_eligible_b <- readr::read_rds(
  here::here("output","data", "data_eligible_b.rds")) 

# read wide vaccine dates data
data_vax_wide <- readr::read_rds(
  here::here("output", "data", "data_wide_vax_dates.rds"))

# read second vax period dates 
second_vax_period_dates <- readr::read_rds(
  here::here("output", "second_vax_period", "data", "second_vax_period_dates.rds")) %>%
  rename(svp_start_date = start_of_period, svp_end_date = end_of_period)

################################################################################
# covariate data
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) %>%
  select(patient_id, subgroup, endoflife_date, midazolam_date, covid_any_date, longres_date)

################################################################################
# apply eligibility criteria in box c ----
data_eligible_c <- data_eligible_b %>%
  left_join(data_vax_wide, 
            by = "patient_id") %>%
  # keep brand of interest 
  # (already applied condition that 1st and 2nd doses are the same)
  mutate(brand = case_when(covid_vax_2_brand %in% "pfizer" ~ "BNT162b2",
                           covid_vax_2_brand %in% "az" ~ "ChAdOx1",
                           TRUE ~ NA_character_)) %>%
  select(-ends_with("_brand")) %>%
  # right join to keep only the jcvi_group:elig_date:region:brands
  # with > n_threshold individuals vaccinated during the period
  right_join(second_vax_period_dates, 
            by = c("jcvi_group", "elig_date", "region")) %>%
  filter(
    # second dose during second vax period
    svp_start_date <= covid_vax_2_date,
    covid_vax_2_date <= svp_end_date) %>%
  select(patient_id, jcvi_group, elig_date, region, ethnicity, 
         covid_vax_2_date, covid_vax_3_date, brand, 
         svp_start_date, svp_end_date) %>%
  droplevels()

# count the number of patients in the extracted data
eligibility_count <- tribble(
  ~description, ~n, ~stage,
  "vax: second dose received during SVP.", n_distinct(data_eligible_c$patient_id), "c-in"
)

################################################################################
# apply eligibility criteria in box d ----

# set seed so that 50:50 split reproducible
set.seed(study_parameters$seed)

data_eligible_d <- data_eligible_a %>%
  # randomly split the unvax 50:50 
  # one group for odd comparisons, one for even
  # so that no overlap in follow-up time across comparisons
  mutate(split = factor(
    rbernoulli(nrow(.), p=0.5), 
    labels = c("odd", "even"))) %>%
  left_join(data_vax_wide %>%
              select(-ends_with("_brand")),
            by = "patient_id") %>%
  left_join(second_vax_period_dates, 
             by = c("jcvi_group", "elig_date", "region")) %>%
  # remove individuals who had received any vaccination before the start of the second vax period
  filter(
    is.na(covid_vax_1_date) | covid_vax_1_date >= svp_start_date
  ) %>%
  select(patient_id, jcvi_group, elig_date, region, ethnicity, 
         covid_vax_1_date, svp_start_date, svp_end_date, split) %>%
  droplevels()

eligibility_count <- eligibility_count %>%
  add_row(
    description = "unvax: unvaccinated at start of SVP",
    n =  n_distinct(data_eligible_d$patient_id),
    stage = "d-in"
  )

################################################################################
# apply eligibility criteria in box e ----

exclusion_e <- function(group) {
  
  # function to be applied in dplyr::filter
  no_evidence_of <- function(cov_date, index_date) {
    is.na(cov_date) | index_date <= cov_date
  }
  
  # define data based on group
  if (group == "vax") {
    data <- data_eligible_c
  } else {
    data <- data_eligible_d
  }
  
  # remove if any covid before start of period
  data <- data %>%
    left_join(data_processed, by = "patient_id") %>%
    filter(
      no_evidence_of(covid_any_date, svp_start_date - weeks(12))) 
  
  eligibility_count_e <- tribble(
    ~description, ~n, 
    glue("{group}: Evidence of COVID before SVP - 12 weeks."), n_distinct(data$patient_id)
  )
  
  # remove if in long-term residential home before start date
  data <- data %>%
    filter(
      no_evidence_of(longres_date, svp_start_date))
  
  eligibility_count_e <- eligibility_count_e %>%
    add_row(
      description = glue("{group}: Evidence of longres before SVP."),
      n =  n_distinct(data$patient_id)
    )
  
  # remove if end of life care before start date
  data <- data %>%
    filter(
      no_evidence_of(endoflife_date, svp_start_date),
      no_evidence_of(midazolam_date, svp_start_date)
    ) %>%
    select(-all_of(names(data_processed)[!names(data_processed) %in% "patient_id"]))
  
  eligibility_count_e <- eligibility_count_e %>%
    add_row(
      description = glue("{group}: Evidence of end of life care before SVP."),
      n =  n_distinct(data$patient_id)
    )
  
  list(data = data, eligibility_count = eligibility_count_e)
    
}

################################################################################

data_eligible_e_vax <- exclusion_e("vax")
data_eligible_e_unvax <- exclusion_e("unvax")

readr::write_rds(
  data_eligible_e_vax[[1]],
  here::here("output", "data", "data_eligible_e_vax.rds"),
  compress = "gz")

readr::write_rds(
  data_eligible_e_unvax[[1]],
  here::here("output", "data", "data_eligible_e_unvax.rds"),
  compress = "gz")

################################################################################
# eligibility count
eligibility_count <- eligibility_count %>%
  bind_rows(
    data_eligible_e_vax[[2]] %>% mutate(stage = "e-ex"), 
    data_eligible_e_unvax[[2]] %>% mutate(stage = "e-ex")
  )

# number of people eligible at each stage ----
eligibility_count_cde <- eligibility_count %>%
  mutate(group = str_extract(description, "\\w+:")) %>%
  arrange(group) %>%
  # round up to nearest 7
  mutate(across(n, ~ceiling_any(.x, to=7))) %>%
  group_by(group) %>%
  mutate(n_removed = lag(n) - n) %>%
  ungroup() %>%
  select(-group)

readr::write_csv(
  eligibility_count_cde,
  here::here("output", "tables", "eligibility_count_cde.csv"))

################################################################################
# for reading into study_definition_tests
data_eligible_e <- bind_rows(
  data_eligible_e_vax[[1]] %>% 
    transmute(patient_id, 
              elig_date, 
              svp_start_date,
              start_1_date = covid_vax_2_date + days(14),
              end_1_date = start_1_date + days(28),
              arm = "vax"),
  data_eligible_e_unvax[[1]] %>% 
    transmute(patient_id, 
              elig_date, 
              svp_start_date,
              start_1_date = svp_start_date + days(14),
              end_1_date = start_1_date + days(56),
              arm = "unvax") 
) %>%
  left_join(data_processed %>% select(patient_id, subgroup),
            by = "patient_id") %>%
  group_by(subgroup) %>%
  mutate(min_elig_date = min(elig_date)) %>%
  ungroup() %>%
  select(-subgroup)

for (k in 2:study_parameters$K) {
  
  data_eligible_e <- data_eligible_e %>%
    mutate(
      !! sym(glue("start_{k}_date")) := !! sym(glue("start_{k-1}_date")) + days(28),
      !! sym(glue("end_{k}_date")) := !! sym(glue("end_{k-1}_date")) + days(28),
    )
  
}

data_eligible_e <- data_eligible_e %>%
  mutate(across(ends_with("_date"), as.POSIXct)) 

readr::write_csv(
  data_eligible_e,
  here::here("output", "data", "data_eligible_e.csv")
)


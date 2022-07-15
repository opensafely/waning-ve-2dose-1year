################################################################################
# process comparisons data
library(tidyverse)
library(lubridate)
library(glue)

################################################################################
## import study_parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# import variable names
model_varlist <- readr::read_rds(
  here::here("analysis", "lib", "model_varlist.rds")
)

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))

################################################################################
# individuals eligible based on box c, d & e criteria 
# arm and split info
data_arm <- bind_rows(
  readr::read_rds(
    here::here("output", "data", "data_eligible_e_vax.rds")) %>%
    rename(arm=brand),
  readr::read_rds(
    here::here("output", "data", "data_eligible_e_unvax.rds")) %>%
    mutate(arm = "unvax")
)  %>%
  select(patient_id, arm, split)

# vars from data_processed
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) %>%
  select(patient_id, subgroup,
         jcvi_group, elig_date, region, 
         dereg_date, death_date,
         starts_with(unname(outcomes)),
         any_of(unname(model_varlist$demographic))) %>%
  # add episode end date for censoring follow-up when the outcome is:
  # covidadmitted, coviddeath, noncoviddeath
  mutate(episode_end_date = postest_date + weeks(12))

# vax data
data_wide_vax_dates <- readRDS(
  here::here("output", "data", "data_wide_vax_dates.rds")) %>%
  select(patient_id, covid_vax_1_date, covid_vax_3_date)

# read data for ever covariates
data_covariates <- arrow::read_feather(
  file = here::here("output", "input_covs.feather")) 

################################################################################
data_all <- data_arm %>%
  # join to covariates data
  left_join(
    data_covariates %>%
      select(patient_id, 
             matches(c("start_\\d+_date", "end_\\d+_date")),
             starts_with("anytest"), asplenia,
             any_of(unname(unlist(model_varlist)))) %>%
      mutate(across(contains("_date"), 
                    ~ floor_date(
                      as.Date(.x, format="%Y-%m-%d"),
                      unit = "days"))) %>%
      # replace end_k_date with end_date_model (end of available hospitalisation data)
      mutate(across(starts_with("end_"), 
                    ~pmin(.x, as.Date(study_parameters$end_date_model)))) 
    , by = "patient_id") %>%
  # join to data_processed
  left_join(
    data_processed, by = "patient_id"
  ) %>%
  # join to vaccines
  left_join(
    data_wide_vax_dates, 
    by = "patient_id"
  ) %>%
  # derive remaining covariates
  mutate(
    
    pregnancy = pregnancy & (sex == "Female") & (age < 50),
    
    immunosuppressed = immunosuppressed | asplenia,
    
    multimorb =
      # as.integer(bmi %in% "Obese III (40+)") +
      as.integer(chd)  +
      as.integer(diabetes) +
      as.integer(cld) +
      as.integer(ckd) +
      as.integer(crd) +
      as.integer(immunosuppressed) +
      as.integer(cns),
    
    multimorb = cut(
      multimorb,
      breaks = c(0, 1, 2, Inf),
      labels=c("0", "1", "2+"),
      right=FALSE)
    
  ) %>%
  mutate(across(test_hist_n,
                ~ factor(case_when(
                  is.na(.x) ~ NA_character_,
                  .x < 1 ~ "0",
                  .x < 2 ~ "1",
                  .x < 3 ~ "2",
                  TRUE ~ "3+"
                )))) %>%
  mutate(subsequent_vax_date = if_else(
    arm %in% "unvax",
    covid_vax_1_date,
    covid_vax_3_date)) %>%
  select(-covid_vax_1_date, -covid_vax_3_date, -asplenia)

readr::write_rds(
  data_all,
  here::here("output", "data", "data_all.rds"),
  compress = "gz"
)

################################################################################
# store min and max fu dates for each subgroup

# create output directory
fs::dir_create(here::here("output", "lib"))

# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

end_K_date <- glue("end_{K}_date")

data_min_max_fu <- data_all %>%
  rename("end_K_date" = end_K_date) %>%
  group_by(subgroup) %>%
  summarise(
    min_fu_date = min(start_1_date),
    max_fu_date = max(end_K_date),
    # round total to nereast 7 for disclosure control
    n = ceiling_any(n(), to=7),
    .groups = "keep"
  ) %>% 
  ungroup() %>%
  mutate(across(max_fu_date,
                ~ pmin(as.Date(study_parameters$end_date), .x)))


# data for release
readr::write_csv(
  data_min_max_fu,
  here::here("output", "lib", "data_min_max_fu.csv")
)


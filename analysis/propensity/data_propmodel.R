################################################################################
#  prepare data for 3rd dose model

# setup ----
library(tidyverse)
library(lubridate)
library(glue)
library(survival)

# source functions
source(here::here("analysis", "functions", "data_process_functions.R"))

# define additional functions
gluec <- function(x) as.character(glue(x))

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# real list of model covariates
model_varlist <- readr::read_rds(
  here::here("analysis", "lib", "model_varlist.rds")
)

# create output directory
outdir <- here::here("output", "propensity", "data")
fs::dir_create(outdir)

# load data ----

# read data from individuals in the 2 dose group who were included in the VE analysis
data_included <- lapply(
  subgroup_labels,
  function(x) {
    comparison <- "both"
    # only one brand used in subgroups 3 and 4
    if (x==3) comparison <- "ChAdOx1"
    if (x==4) comparison <- "BNT162b2"
    readr::read_rds(here::here("output", "tte", "data", glue("data_tte_{comparison}_{x}_coviddeath.rds"))) %>%
      filter(arm != "unvax") %>%
      distinct(patient_id)
  }
) %>%
  bind_rows()

# read data with baseline covariates and outcomes
data_all <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) %>%
  # only keep those in both datasets
  inner_join(data_included, by = "patient_id")

# covariates
data_covs <- data_all %>%
  select(
    patient_id, 
    subgroup, jcvi_group, elig_date, region,
    vax2_date = start_1_date,
    vax2_brand = arm, 
    all_of(unname(unlist(model_varlist)))
  ) %>%
  mutate(across(vax2_brand, ~factor(.x)))

# time to event
data_tte <- data_all %>%
  transmute(
    patient_id, 
    end_fup_date = pmin(dereg_date, death_date, end_12_date, subsequent_vax_date, na.rm=TRUE),
    day = as.integer(end_fup_date - start_1_date),
    status = !is.na(subsequent_vax_date) & subsequent_vax_date <= end_fup_date
  ) 

# load extracted data
data_extract <- arrow::read_feather(here::here("output", "input_prop.feather")) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(
    across(
      c(contains("_date")), 
      ~ floor_date(as.Date(., format="%Y-%m-%d"), unit = "days")
    )
  ) 


# define recurring variables in extracted data
recurring_vars <- c(
  "postest", 
  "admitted_planned", "admitted_unplanned", 
  # covidunplanned must come after unplanned here
  "admitted_covidunplanned"
)

# sort out dummy data ----

## if running on dummy data make it more sensible
# i.e. enforce the following:
# - all dates after start_1_date
# - recurring_var_2_date > recurring_var_1_date
# - discharged_var_1_date > admitted_var_1_date (but allow for some errors here as in real data)
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  
  data_extract <- local({
    
    # need to redefine gluec() inside local({}) otherwise get errors
    gluec <- function(x) as.character(glue(x))
    extract_tmp <- data_extract
    
    for (v in recurring_vars) {
      
      indices <- extract_tmp %>% select(starts_with(v)) %>% names() %>% 
        str_extract(., "\\d+") %>% as.integer() %>% sort()
      
      # loop over the number of times each variable recurs
      for (i in indices) {
        
        if (i == 0) {
          # ensure all dates before start_1_date
          extract_tmp <- extract_tmp %>%
            mutate(
              across(
                matches(gluec("{v}_{i}_date")), 
                ~if_else(
                  .x < start_1_date,
                  as.Date(.x),
                  as.Date(NA_character_)
                )))
        } 
        
        if (i == 1) {
          if (str_detect(v, "admitted")) shift <- 0 else shift <- -28
          # ensure all dates are after start_1_date - shift
          extract_tmp <- extract_tmp %>%
            mutate(
              across(
                matches(gluec("{v}_{i}_date")),
                ~if_else(
                  start_1_date + shift < .x,
                  as.Date(.x),
                  as.Date(NA_character_)
                )))
        }
        
        if (i > 1) {
          # define name of the variable before the current one
          v_before <- sym(gluec("{v}_{i-1}_date"))
          # make the current variable missing if it is not after the one before
          extract_tmp <- extract_tmp %>%
            mutate(
              across(
                matches(gluec("{v}_{i}_date")),
                ~if_else(
                  !! v_before < .x,
                  as.Date(.x),
                  as.Date(NA_character_)
                )))
        }
        
        if (str_detect(v, "admitted") & (v != "admitted_covidunplanned")) {
          # make sure discharges occur after admissions
          v_admitted <- gluec("{v}_{i}_date")
          v_discharged <- str_replace(v_admitted, "admitted", "discharged")
          extract_tmp <- extract_tmp %>%
            # should be 1:60, but use -1:60 to introduce some errors
            mutate(!! sym(v_discharged) := !! sym(v_admitted) + runif(nrow(.), -1, 60))
        }
        
        if (v == "admitted_covidunplanned") {
          extract_tmp <- extract_tmp %>%
            mutate(
              across(
                gluec("admitted_covidunplanned_{i}_date"),
                ~if_else(
                  is.na(.x),
                  as.Date(NA_character_),
                  !!sym(gluec("admitted_unplanned_{i}_date"))
                )
              )
            ) %>%
            mutate(
              across(
                gluec("discharged_covidunplanned_{i}_date"),
                ~if_else(
                  is.na(!!sym(gluec("admitted_covidunplanned_{i}_date"))),
                  as.Date(NA_character_),
                  !!sym(gluec("discharged_unplanned_{i}_date"))
                )
              )
            )
        }
        
      }
      
    }
    
    return(extract_tmp)
    
  })
  
}

#### process data ----
  
data_processed <- data_extract %>%
  # only keep those in both datasets
  inner_join(data_included, by = "patient_id") %>%
  # join end_fup_date
  left_join(data_tte %>% select(patient_id, end_fup_date), by = "patient_id") %>%
  # remove dates after end_fup_date as not needed
  mutate(
    across(
      ends_with("_date"), 
      ~if_else(
        .x <= end_fup_date,
        .x,
        as.Date(NA_character_)
      )
    )
  )

# clean up
rm(data_all, data_included, data_extract)

# define static vars as used frequently in following chunks
static_vars <- c("patient_id", "start_1_date")

## endoflife:
data_endoflife <- data_processed %>%
  select(all_of(static_vars), endoflife_date) %>%
  transmute(
    patient_id, 
    endoflife = !is.na(endoflife_date),
    day = if_else(
      endoflife,
      as.integer(endoflife_date - start_1_date),
      0L
    )
  )

## postest:
data_postest <- data_processed %>%
  select(all_of(static_vars), starts_with("postest")) %>%
  pivot_longer(cols = -all_of(static_vars), values_drop_na = TRUE) %>%
  select(all_of(static_vars), postest_start_date = value) %>%
  # end date for a "positive test episode" is 30 days after date of positive test
  mutate(postest_end_date = postest_start_date + 30) %>%
  # get the date of the next positive test
  group_by(patient_id) %>%
  mutate(lead_postest_start_date = lead(postest_start_date)) %>%
  ungroup() %>%
  # if the next positive test is less than 30 days after the current one, remove the end date
  mutate(
    across(
      postest_end_date, 
      ~if_else(
        !is.na(lead_postest_start_date) & (lead_postest_start_date <= .x),
        as.Date(NA_character_),
        .x
      ))) %>%
  select(-lead_postest_start_date) %>%
  # transform data to long format
  pivot_longer(cols = -all_of(static_vars), values_drop_na = TRUE) %>%
  # postest = TRUE for start dates, = FALSE for end dates
  mutate(postest = name == "postest_start_date") %>%
  # arrange events (start and end dates) by their date, within individuals
  group_by(patient_id) %>%
  arrange(value, .by_group = TRUE) %>%
  # calculate a "moving sum" of the current row and the previous row, 
  # to flag instances where there are two start dates in a row
  mutate(cumulsum = stats::filter(postest, rep(1,2), sides = 1)) %>%
  ungroup() %>% 
  # remove start dates if preceding row is also a start date (within patients)
  filter(is.na(cumulsum) | cumulsum == 1) %>%
  transmute(
    patient_id,
    day = as.integer(value - start_1_date),
    postest
  )
  
## clean hospitalisation data
# Need to sort the following issues:
# issue1: admission_*_x_date > discharge_*_x_date: 
#         stay of length 1 starting on admission date
# issue2: admission_*_x_date == discharge_*_x_date: 
#         stay of length 1
# issue3: admission_*_x_date <= discharge_*_{x-1}_date: 
#         stay continuous lasting from admission_*_x_date until discharge_*_{x+1}_date

hosp_vars <- recurring_vars[str_detect(recurring_vars, "admitted")]
issues <- list()
data_hosp_processed <- data_processed %>% 
  select(all_of(static_vars), starts_with(c("admitted", "discharged")))

# summarise and fix the above issues
for (v in hosp_vars) {
  
  indices <- data_hosp_processed %>% select(starts_with(v)) %>% names() %>% 
    str_extract(., "\\d+") %>% as.integer() %>% sort()
  
  v_type <- str_remove(v, "admitted_")
  v_discharged <- str_replace(v, "admitted", "discharged")
  
  issues[[v]] <- list(issue1 = integer(), issue2 = integer(), issue3 = integer())
  
  # summarise the number of times each issue occurs
  for (i in indices) {
    
    # issue1
    issues[[v]][["issue1"]] <- sum(
      issues[[v]][["issue1"]],
      data_hosp_processed %>%
        filter(!!sym(gluec("{v}_{i}_date")) > !!sym(gluec("{v_discharged}_{i}_date"))) %>%
        nrow()
    )
    
    # issue2
    issues[[v]][["issue2"]] <- sum(
      issues[[v]][["issue2"]],
      data_hosp_processed %>%
      filter(!!sym(gluec("{v}_{i}_date")) == !!sym(gluec("{v_discharged}_{i}_date"))) %>%
      nrow()
    )
    
    # issue3
    if (i>0) {
      issues[[v]][["issue3"]] <- sum(
        issues[[v]][["issue3"]],
      data_hosp_processed %>%
        filter(!!sym(gluec("{v}_{i}_date")) <= !!sym(gluec("{v_discharged}_{i-1}_date"))) %>%
        nrow()
      )
    }
    
  }
  
  # fix issues 1&2
  # note: this isn't a great fix for issue 1, but I can't think of any better 
  # way without excluding patients with this issue, which we don't want to do
  for (i in indices) {
    data_hosp_processed <- data_hosp_processed %>%
      mutate(
        across(
          gluec("{v_discharged}_{i}_date"),
          ~if_else(
            # discharged before admitted
            .x <= !!sym(gluec("{v}_{i}_date")),
            !!sym(gluec("{v}_{i}_date")),
            .x
            )
          )
        ) 
  }
  
  # add 1 day to each discharge date so that:
  # - someone admitted and discharged on the same day has a stay of length 1
  # - someone someone admitted the day after discharge from a previous
  #   hospitalisation has a continuous stay across both hospitalisations
  # note: do this before sorting issue 3, otherwise could reintroduce issue 3
  data_hosp_processed <- data_hosp_processed %>%
    mutate(across(starts_with(v_discharged), ~.x+1))
  
  # fix issue 3
  for (i in indices[-1]) {
    
    data_hosp_processed <- data_hosp_processed %>%
      mutate(
        across(
          c(gluec("{v_discharged}_{i-1}_date"), gluec("{v}_{i}_date")),
          ~ if_else(
            !!sym(gluec("{v}_{i}_date")) <= !! sym(gluec("{v_discharged}_{i-1}_date")),
            as.Date(NA_character_),
            .x
          )
          )
        )
      
  }
  
}

# print the issues summary
print(lapply(issues, function(x) unlist(x)))

# derive data_hosp:
data_hosp <- list()
for (v in hosp_vars) {
  
  v_type <- str_remove(v, "admitted_")
  data_hosp[[v]] <- data_hosp_processed %>%
    select(all_of(static_vars), matches(glue("\\w+_{v_type}_\\d+_date"))) %>%
    pivot_longer(
      cols = -all_of(static_vars),
      values_drop_na = TRUE
    ) %>%
    transmute(
      patient_id,
      day = as.integer(value - start_1_date),
      # inhosp_* = TRUE if an admitted date, FALSE if a discharged date
      !!sym(gluec("inhosp_{v_type}")) := str_detect(name, "^admitted")
    )
  
}

# join datasets of time dependent covariates
by_vars <- c("patient_id", "day")
data_tdc <- data_endoflife %>%
  full_join(data_postest, by = by_vars) %>%
  full_join(data_hosp$admitted_planned, by = by_vars) %>%
  full_join(data_hosp$admitted_unplanned, by = by_vars) %>%
  full_join(data_hosp$admitted_covidunplanned, by = by_vars) 

# clean up
rm(data_postest, data_hosp, data_endoflife, data_hosp_processed)

cat("check for days with >1 rows:\n")
data_tdc %>% 
  group_by(patient_id, day) %>%
  mutate(n=n()) %>% 
  filter(n>1) %>%
  nrow()

# fill missing values from joins
data_filled <- data_tdc %>%
  select(patient_id, day, everything()) %>%
  arrange(patient_id, day) %>%
  # create `sequence` of dates within patients
  group_by(patient_id) %>%
  mutate(sequence = row_number()) %>%
  ungroup() %>%
  # if vairables are missing when sequence = 1, fill with FALSE
  # (if it were TRUE they would have been excluded previously)
  mutate(
    across(
      -c(patient_id, day), 
      ~if_else((sequence == 1) & (is.na(.x)), FALSE, as.logical(.x))
      )
    ) %>%
  select(-sequence) %>%
  # within each patient, fill missing values with the last non missing value
  group_by(patient_id) %>%
  fill(-c(patient_id, day), .direction = "down") %>%
  ungroup()

# for patients whose first entry is after day 0, impute day 0
# we know all time dependent covariates must be FALSE, otherwise they would 
# have an entry or have been excluded already
data_day0 <- data_filled %>%
  anti_join(
    data_filled %>% filter(day <= 0) %>% distinct(patient_id), 
    by = "patient_id"
    ) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  mutate(day = 0) %>%
  mutate(across(-c(patient_id, day), ~FALSE)) 

# bind datasets
data_final <- bind_rows(
  data_day0,
  data_filled
) %>%
  arrange(patient_id, day) 

# clean up
rm(data_tdc, data_day0, data_filled)

# derive final data for cox model
data_propmodel <- data_tte %>% 
  select(-end_fup_date) %>%
  # derive tstart and tstop
  tmerge(
    data1 = .,
    data2 = .,
    id = patient_id,
    tstart = 0,
    tstop = day,
    ind_outcome = event(day, status)
  ) %>% 
  # derive time dependent covariates
  tmerge(
    data1 = .,
    data2 = data_final,
    id = patient_id,
    # time dependent covariates
    endoflife = tdc(day, endoflife),
    postest = tdc(day, postest),
    inhosp_planned = tdc(day, inhosp_planned),
    inhosp_unplanned = tdc(day, inhosp_unplanned),
    inhosp_covidunplanned = tdc(day, inhosp_covidunplanned)
  ) %>%
  as_tibble() %>%
  select(-c(day, status)) %>%
  # join baseline covariates
  left_join(data_covs, by = "patient_id")

# save dataset
write_rds(
  data_propmodel,
  file.path(outdir, "data_propmodel.rds"),
  compress = "gz"
)

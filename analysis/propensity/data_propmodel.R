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

# time to event
data_tte <- readr::read_rds(here::here("output", "data", "data_all.rds")) %>%
  # only keep those in both datasets
  inner_join(data_included, by = "patient_id") %>%
  transmute(
    patient_id, 
    end_fup_date = pmin(dereg_date, death_date, end_12_date, subsequent_vax_date, na.rm=TRUE),
    day = as.integer(end_fup_date - start_1_date),
    status = as.integer(!is.na(subsequent_vax_date) & subsequent_vax_date <= end_fup_date)
  ) 

# load extracted data
data_extract <- arrow::read_feather(here::here("output", "input_prop.feather")) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(
    across(
      c(contains("_date")), 
      ~ floor_date(as.Date(., format="%Y-%m-%d"), unit = "days")
    )
  ) %>%
  # only keep those in both datasets
  inner_join(data_included, by = "patient_id") %>%
  # for now only flag a cancer diagnosis in the 3 years before start_1_date
  mutate(cancer_0_date = max(cancerhosp_0_date, cancerprimarycare_0_date)) %>%
  mutate(
    across(
      cancer_0_date, 
      ~if_else(
        cancer_0_date < start_1_date - 3*365,
        as.Date(NA_character_), 
        .x
        )
      )
    ) %>%
  select(-matches("cancer\\w+_\\d+_date"))

# process_data ----

data_long <- data_extract  %>%
  pivot_longer(
    cols = -c(patient_id, start_1_date),
    values_drop_na = TRUE
    ) %>%
  mutate(across(name, ~str_remove(.x, "_\\d+"))) %>%
  mutate(across(name, ~str_remove(.x, "_date$"))) %>%
  mutate(day = as.integer(value - (start_1_date-1))) %>%
  select(-c(start_1_date, value)) 

# add postest stop times
data_long <- data_long %>%
  mutate(across(name, ~str_replace(.x, "postest", "postest_start"))) %>%
  bind_rows(
    data_long %>% 
      filter(name == "postest") %>%
      mutate(
        name = "postest_end",
        day = day+30
        )
  ) %>%
  arrange(patient_id, day)

data_timevarying <- data_tte %>%
  select(patient_id) %>%
  arrange(patient_id) %>%
  tmerge(
    # initialise tmerge run
    data1 = .,
    data2 = data_tte, 
    id = patient_id,
    tstart = 0,
    tstop = day,
    ind_outcome = event(day, status)
  ) %>% 
  tmerge(
    data1 = .,
    data2 = data_long,
    id = patient_id,
    status_postest = tdc(
      day, 
      case_when(
        name=="postest_start" ~ 1L, 
        name=="postest_end" ~ 0L, 
        TRUE ~NA_integer_
        )
      ),
    status_endoflife = event(day, as.integer(name == "endoflife")),
    status_cancer = event(day, as.integer(name == "cancer")),
    status_planned = tdc(
      day, 
      case_when(
        name=="admitted_planned" ~ 1L, 
        name=="discharged_planned" ~ 0L, 
        TRUE ~NA_integer_
      )
    ),
    status_unplanned = tdc(
      day, 
      case_when(
        name=="admitted_unplanned" ~ 1L, 
        name=="discharged_unplanned" ~ 0L, 
        TRUE ~NA_integer_
      )
    ),
    status_covidunplanned = tdc(
      day, 
      case_when(
        name=="admitted_covidunplanned" ~ 1L, 
        name=="discharged_covidunplanned" ~ 0L, 
        TRUE ~NA_integer_
      )
    )
  ) %>% 
  # if it's missing it must be absent, otherwise they would have been excluded
  mutate(across(starts_with("status_"), ~replace_na(data = .x, replace = 0L))) %>%
  as_tibble()

# save dataset
write_rds(
  data_timevarying,
  file.path(outdir, "data_timevarying.rds"),
  compress = "gz"
)

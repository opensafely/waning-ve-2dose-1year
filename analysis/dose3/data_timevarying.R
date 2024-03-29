################################################################################
#  prepare data for 3rd dose model

# setup ----
library(tidyverse)
library(lubridate)
library(glue)
library(survival)

# source functions
source(here::here("analysis", "functions", "data_process_functions.R"))
source(here::here("analysis", "functions", "data_properties.R"))

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
outdir <- here::here("output", "dose3", "data")
fs::dir_create(outdir)

# load data ----

# read data from individuals in the 2 dose group who were included in the VE analysis
data_included <- lapply(
  subgroup_labels,
  function(x) {
    df <- list()
    # only one brand used in subgroups 3 and 4
    if (x %in% c(1, 2, 3)) {
      df[[1]] <- readr::read_rds(here::here("output", "tte", "data", glue("data_tte_ChAdOx1_{x}_coviddeath.rds")))
    } 
    if (x %in% c(1, 2, 4)) {
      df[[2]] <- readr::read_rds(here::here("output", "tte", "data", glue("data_tte_BNT162b2_{x}_coviddeath.rds")))
    } 
     bind_rows(df) %>%
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
data_extract0 <- arrow::read_feather(here::here("output", "input_timevarying.feather")) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(
    across(
      c(contains("_date")), 
      ~ floor_date(as.Date(., format="%Y-%m-%d"), unit = "days")
    )
  )

# save data properties
data_properties(
  data = data_extract0,
  path = outdir
)

data_extract <- data_extract0 %>%
  # only keep those in both datasets
  inner_join(data_included, by = "patient_id") %>%
  # for now only flag a cancer diagnosis in the 3 years before start_1_date
  mutate(cancer_0_date = pmax(cancerhosp_0_date, cancerprimarycare_0_date, na.rm=TRUE)) %>%
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

# save data properties
data_properties(
  data = data_extract,
  path = outdir
)


# process_data ----

data_long <- data_extract  %>%
  pivot_longer(
    cols = -c(patient_id, start_1_date),
    values_drop_na = TRUE
    ) %>%
  # remove number and _date from names
  mutate(across(name, ~str_remove(.x, "_\\d+"))) %>%
  mutate(across(name, ~str_remove(.x, "_date$"))) %>%
  # rescale to days since start date (days=1 on start date)
  mutate(day = as.integer(value - (start_1_date-1))) %>%
  select(-c(start_1_date, value)) 

# define postest_end as 30 days after each postest
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
  bind_rows(
    # adding this dummy category is a hack to make sure these is always
    # a period for every patient starting at start=0
    data_tte %>% distinct(patient_id) %>% mutate(name = "dummy", day = 0)
    ) %>%
  arrange(patient_id, day)

# save data properties
data_properties(
  data = data_long,
  path = outdir
)

data_timevarying <- data_tte %>%
  select(patient_id) %>%
  arrange(patient_id) %>%
  tmerge(
    # initialise tmerge run
    data1 = .,
    data2 = data_tte, 
    id = patient_id,
    # because the earliest date will be for cancer, 
    # which must be <=3*365 days before start_1_date 
    tstart = - 3*365 - 1,
    tstop = day,
    ind_outcome = event(day, status)
  ) %>% 
  # define time-varying variables
  tmerge(
    data1 = .,
    data2 = data_long,
    id = patient_id,
    dummy = event(day, name == "dummy"),
    status_endoflife = cumtdc(day, as.integer(name == "endoflife")),
    status_cancer = cumtdc(day, as.integer(name == "cancer")),
    status_postest = tdc(
      day, 
      case_when(
        name=="postest_start" ~ 1L, 
        name=="postest_end" ~ 0L, 
        TRUE ~NA_integer_
      )
    ),
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
    ),
    # initialise all with 0, otherwise they would have been excluded
    options=list(tdcstart=0L)
  ) %>%
  select(-dummy) %>%
  # get rid of all periods that start before start_1_date
  filter(tstart >= 0) %>%
  as_tibble()

# save data properties
data_properties(
  data = data_timevarying,
  path = outdir
)

# save dataset
write_rds(
  data_timevarying,
  file.path(outdir, "data_timevarying.rds"),
  compress = "gz"
)

######################################

# This script:
# - creates dummy data for study_definition.py

######################################

library(tidyverse)
library(lubridate)
library(glue)

source(here::here("analysis", "functions", "dummy_data_functions.R"))

set.seed(5476)

# date vars 
# set these to have occured since start of pandemic
date_vars_recent <- c("positive_test_0_date", 
                      "primary_care_covid_case_0_date", 
                      "covidadmitted_0_date",
                      "covidemergency_0_date",
                      "death_date",
                      "longres_date",
                      "endoflife_date", 
                      "midazolam_date",
                      "coviddeath_date", 
                      "dereg_date")

jcvi_group_patterns <- readr::read_csv(here::here("analysis", "lib", "jcvi_groups.csv")) %>%
  mutate(across(definition, ~str_extract(.x, "age_. >=\\d{2}"))) %>%
  # add dummy conditions for groups 1 and 6, as longres and atrisk data not available here (done correctly in real data)
  mutate(across(definition, ~case_when(group=="01" ~ "age_1 >=85", 
                                       group=="04b" ~ "age_1 >=68",
                                       group=="06" ~ "age_1 >=62",
                                       !is.na(.x) ~ .x,
                                       TRUE ~ "TRUE")))

# conditions for eligibility dates
elig_date_patterns <- readr::read_csv(here::here("analysis", "lib", "elig_dates.csv")) %>%
  mutate(across(description, ~str_replace_all(.x, "p=", "p=="))) %>%
  mutate(across(description, ~str_replace_all(.x, "OR", "|"))) %>%
  mutate(across(description, ~str_replace_all(.x, "AND", "&"))) %>%
  mutate(across(description, ~str_replace_all(.x, "DEFAULT", "TRUE")))

dummy_data_elig <- tibble(patient_id = 1:n) %>%
  mutate(age_1 = as.integer(runif(nrow(.), 16, 90), 0),
         age_2 = age_1,
         imd = sample(
           x = seq(100L,32100L,100L),
           size = n,
           replace = TRUE)) %>%
  var_category(sex, categories = c("F", "M")) %>%
  var_category(ethnicity_6, 
               categories = c(as.character(1:5), NA_character_),
               ratios = c(rep(0.99/5,5), 0.01)) %>%
  var_category(ethnicity_6_sus, 
               categories = c(as.character(1:5), NA_character_),
               ratios = c(rep(0.99/5,5), 0.01)) %>%
  var_category(region, categories = regions$region, ratios = regions$ratio) %>%
  # binary vars for exclusion criteria
  mutate(hscworker = rbernoulli(n = nrow(.), p=0.01)) %>%
  mutate(housebound = rbernoulli(n = nrow(.), p=0.01)) %>%
  # jcvi_group
  var_category(
    name = jcvi_group, 
    categories = jcvi_group_patterns$group, 
    conditions = jcvi_group_patterns$definition
    ) %>%
  # elig_date
  var_category(
    name = elig_date,
    categories = elig_date_patterns$date,
    conditions = elig_date_patterns$description) 

# fix vaccine dates so that they have roughly correct distribution
# fix vaccine dates so that they have roughly correct distribution
dummy_data_vax <- dummy_data_elig %>%
  mutate(
    covid_vax_pfizer_1_date = as.Date(elig_date) + days(round(rnorm(nrow(.), mean = 10, sd = 3))),
    covid_vax_az_1_date = as.Date(elig_date) + days(round(rnorm(nrow(.), mean = 10, sd = 3))),
    covid_vax_moderna_1_date = as.Date(elig_date) + days(round(rnorm(nrow(.), mean = 10, sd = 3))),
    covid_vax_disease_1_date = as.Date(elig_date) + days(round(rnorm(nrow(.), mean = 10, sd = 3)))) %>%
  mutate(
    vaccine_1_type = sample(
      x = c("pfizer", "az", "moderna", "disease", "none"), 
      size = nrow(.),
      replace = TRUE,
      prob = c(0.4, 0.4, 0.05, 0.05, 0.1)
    ),
    missing_pfizer_2 = rbernoulli(nrow(.), p=0.05),
    missing_az_2 = rbernoulli(nrow(.), p=0.05),
    missing_moderna_2 = rbernoulli(nrow(.), p=0.05),
    missing_pfizer_3 = rbernoulli(nrow(.), p=0.9),
    missing_az_3 = rbernoulli(nrow(.), p=0.9),
    missing_moderna_3 = rbernoulli(nrow(.), p=0.9)
  ) %>%
  mutate(across(covid_vax_pfizer_1_date, 
                ~if_else(
                  vaccine_1_type %in% "pfizer",
                  .x,
                  NA_Date_))) %>%
  mutate(across(covid_vax_az_1_date, 
                ~if_else(
                  vaccine_1_type %in% "az",
                  .x,
                  NA_Date_))) %>%
  mutate(across(covid_vax_moderna_1_date, 
                ~if_else(
                  vaccine_1_type %in% "moderna",
                  .x,
                  NA_Date_))) %>%
  mutate(across(covid_vax_disease_1_date, 
                ~if_else(
                  vaccine_1_type %in% "disease",
                  .x,
                  NA_Date_))) %>%
  mutate(across(matches("covid_vax_\\w+_1_date"),
                ~ if_else(
                  vaccine_1_type %in% "none",
                  NA_Date_,
                  .x
                ))) %>%
  mutate(
    covid_vax_pfizer_2_date = covid_vax_pfizer_1_date + days(round(rnorm(nrow(.), mean = 10*7, sd = 3))),
    covid_vax_az_2_date = covid_vax_az_1_date + days(round(rnorm(nrow(.), mean = 10*7, sd = 3))),
    covid_vax_moderna_2_date = covid_vax_moderna_1_date + days(round(rnorm(nrow(.), mean = 10*7, sd = 3))),
  ) %>%
  mutate(across(covid_vax_pfizer_2_date, 
                ~if_else(
                  missing_pfizer_2,
                  NA_Date_,
                  .x))) %>%
  mutate(across(covid_vax_az_2_date, 
                ~if_else(
                  missing_az_2,
                  NA_Date_,
                  .x))) %>%
  mutate(across(covid_vax_moderna_2_date, 
                ~if_else(
                  missing_moderna_2,
                  NA_Date_,
                  .x))) %>%
  mutate(
    covid_vax_pfizer_3_date = covid_vax_pfizer_2_date + days(round(rnorm(nrow(.), mean = 6*4*7, sd = 7))),
    covid_vax_az_3_date = covid_vax_az_2_date + days(round(rnorm(nrow(.), mean = 6*4*7, sd = 7))),
    covid_vax_moderna_3_date = covid_vax_moderna_2_date + days(round(rnorm(nrow(.), mean = 6*4*7, sd = 7))),
  ) %>%
  select(-starts_with("missing"), -vaccine_1_type) 

### from dummy_data_covs
dummy_data_covs <- dummy_data_vax %>%
  # date vars recent
  bind_cols(
    pmap(
      list(a = date_vars_recent, 
           b = rep(0.2, length(date_vars_recent))),
      function(a,b) 
        var_date(
          .data = ., 
          name = !! a,
          incidence = b,
          earliest=study_parameters$pandemic_start,
          latest=study_parameters$end_date,
          keep_vars = FALSE
        ))) %>%
  # add death_date if coviddeath_date
  mutate(across(death_date, 
                ~if_else(
                  !is.na(coviddeath_date), 
                  coviddeath_date,
                  .x))) %>%
  mutate(across(contains("_date"), as.POSIXct)) %>%
  mutate(across(ends_with("date"), as.POSIXct)) %>%
  mutate(across(c(ethnicity_6, ethnicity_6_sus, jcvi_group, region, sex),
                as.factor)) %>%
  droplevels()

# final dummy data
dummy_data <- dummy_data_covs

arrow::write_feather(dummy_data, here::here("analysis", "dummy_data_vax.feather"))

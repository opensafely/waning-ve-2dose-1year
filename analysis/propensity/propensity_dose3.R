################################################################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

## source functions
source(here::here("analysis", "functions", "data_process_functions.R"))

## load extracted data
extract <- arrow::read_feather(here::here("output", "input_prop.feather")) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(across(c(contains("_date")), 
                ~ floor_date(
                  as.Date(., format="%Y-%m-%d"),
                  unit = "days")))

recurring_vars <- c(
  "postest", "admitted_planned", "admitted_unplanned", 
  # covidunplanned must come after unplanned here
  "admitted_covidunplanned"
  )

## if running on dummy data make it more sensible
# i.e. enforce the following:
# - all dates after start_1_date
# - recurring_var_2_date > recurring_var_1_date
# - discharged_var_1_date > admitted_var_1_date
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  
  for (v in recurring_vars) {
    indices <- extract %>% select(starts_with(v)) %>% names() %>% 
      str_extract(., "\\d+") %>% as.integer() %>% sort()
    # loop over the number of times each variable recurs
    for (i in indices) {
      
      if (i == 0) {
        # ensure all dates before start_1_date
        extract <- extract %>%
          mutate(
            across(
              matches(glue("{v}_{i}_date")), 
              ~if_else(
                .x < start_1_date,
                as.Date(.x),
                as.Date(NA_character_)
              )))
      } 
      
      if (i == 1) {
        if (str_detect(v, "admitted")) shift <- 0 else shift <- -28
        # ensure all dates are after start_1_date - shift
        extract <- extract %>%
          mutate(
            across(
              matches(glue("{v}_{i}_date")), 
              ~if_else(
                start_1_date + shift < .x,
                as.Date(.x),
                as.Date(NA_character_)
              )))
      } 
      
      if (i > 1) {
        # define name of the variable before the current one
        v_before <- sym(as.character(glue("{v}_{i-1}_date")))
        # make the current variable missing if it is not after the one before
        extract <- extract %>%
          mutate(
            across(
              matches(glue("{v}_{i}_date")), 
              ~if_else(
                !! v_before < .x,
                as.Date(.x),
                as.Date(NA_character_)
              )))
      }
      
      if (str_detect(v, "admitted") & (v != "admitted_covidunplanned")) {
        # make sure discharges occur after admissions
        v_admitted <- as.character(glue("{v}_{i}_date"))
        v_discharged <- str_replace(v_admitted, "admitted", "discharged")
        extract <- extract %>%
          mutate(!! sym(v_discharged) := !! sym(v_admitted) + runif(nrow(.), 1, 60)) 
      } 
      
      if (v == "admitted_covidunplanned") {
        extract <- extract %>%
          mutate(
            across(
              glue("admitted_covidunplanned_{i}_date"),
              ~if_else(
                is.na(.x),
                as.Date(NA_character_),
                !!sym(as.character(glue("admitted_unplanned_{i}_date")))
              )
            )
          ) %>%
          mutate(
            across(
              glue("discharged_covidunplanned_{i}_date"),
              ~if_else(
                is.na(!!sym(as.character(glue("admitted_covidunplanned_{i}_date")))),
                as.Date(NA_character_),
                !!sym(as.character(glue("discharged_unplanned_{i}_date")))
              )
            )
          )
      }
      
    }
  }
  
}

# status on start_1_date
# loop over hospitalisations:
data_start_1_date <- extract
for (v in recurring_vars[str_detect(recurring_vars, "admitted")]) {
  
  v_type <- str_remove(v, "admitted_")
  v_discharged <- as.character(glue("discharged_{v_type}_0_date"))
  data_start_1_date <- data_start_1_date %>%
    mutate(
      !!sym(as.character(glue("inhosp_{v_type}"))) :=
        !is.na(!!sym(as.character(glue("{v}_0_date")))) &
        !is.na(!!sym(v_discharged)) &
        # still counted as in hospital on day of discharge
        start_1_date <= !!sym(v_discharged)   
    )
  
}

data_start_1_date <- data_start_1_date %>%
  transmute(
    patient_id, day = 0, 
    inhosp_planned, inhosp_unplanned, inhosp_covidunplanned,
    # eol care initiated
    eol = !is.na(endoflife_date) &
      (endoflife_date < start_1_date),
    # positive test in past 30 days
    recentpostest = !is.na(postest_1_date) &
      (postest_1_date < start_1_date) & 
      as.integer(start_1_date - postest_1_date) < 30
  )

## each time status updates after start_1_date

# transform data to long format
static_vars <- c("patient_id", "start_1_date")
data_postest_long <- extract %>%
  select(all_of(static_vars), starts_with("postest")) %>%
  pivot_longer(cols = -all_of(static_vars), values_drop_na = TRUE) 

data_postest <- data_postest_long %>%
  select(all_of(static_vars), postest_start_date = value) %>%
  # end date for a "positive test episode" is 30 days after date of positive test
  mutate(postest_end_date = postest_start_date+30) %>%
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
  # postest=1 for start dates, =0 for end dates
  mutate(postest = name == "postest_start_date") %>%
  # arrange events (start and end dates) by their date, within individuals
  group_by(patient_id) %>%
  arrange(value, .by_group = TRUE) %>%
  # calculate a "moving sum" of the current row and the previous row, 
  # to flag instances where there are two start dates in a row
  mutate(cumulsum = stats::filter(postest, rep(1,2), sides = 1)) %>%
  ungroup() %>% 
  # remove start dates that are not separated by end dates
  filter(is.na(cumulsum) | cumulsum == 1) 

# check that there's never a new positive test less than 30 days after the last one
data_postest %>%
  filter(postest) %>%
  group_by(patient_id) %>%
  mutate(timesincelast_postest = c(NA_integer_, diff(value))) %>%
  ungroup() %>%
  filter(timesincelast_postest < 30)

  
  

bind_rows(
  dat_postest_long %>% mutate(name = "postest_date"),
  dat_postest_long %>% mutate(name = "postest_plus30_date", value = value+30)
  ) %>%
  arrange(patient_id, value)

data_updated <- bind_rows(
  extract %>%
    select(all_of(static_vars), starts_with("postest")) %>%
    pivot_longer(cols = -all_of(static_vars), values_drop_na = TRUE) %>%
    mutate(name = "postest_date") %>%
    mutate(value_30 = value+30) %>%
    pivot_longer
    
    pivot_longer(
      cols = -c(patient_id, start_1_date),
      values_drop_na = TRUE
    )
)
  
  
  

## process recurring variables relative to start_1_date

extract %>%
  select(patient_id, start_1_date, matches("\\w+_planned_\\d_date")) %>%
  pivot_longer(cols = -c(patient_id), values_drop_na = TRUE) %>%
  mutate(across(value, ~if_else(str_detect(name, "^discharged"), .x+1, .x))) %>%
  mutate(inhosp_unplanned = case_when(
    name == start_1_date ~ 0L,
    
  ))
  


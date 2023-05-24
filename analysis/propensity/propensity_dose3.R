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

gluec <- function(x) as.character(glue(x))

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
        extract <- extract %>%
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
        extract <- extract %>%
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
        extract <- extract %>%
          # should be 1:60, but use -1:60 to introduce some errors
          mutate(!! sym(v_discharged) := !! sym(v_admitted) + runif(nrow(.), -1, 60)) 
      } 
      
      if (v == "admitted_covidunplanned") {
        extract <- extract %>%
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
  
}

# # status on start_1_date
# # loop over hospitalisations:
# data_start_1_date <- extract
# for (v in recurring_vars[str_detect(recurring_vars, "admitted")]) {
#   
#   v_type <- str_remove(v, "admitted_")
#   v_discharged <- as.character(glue("discharged_{v_type}_0_date"))
#   data_start_1_date <- data_start_1_date %>%
#     mutate(
#       !!sym(as.character(glue("inhosp_{v_type}"))) :=
#         !is.na(!!sym(as.character(glue("{v}_0_date")))) &
#         !is.na(!!sym(v_discharged)) &
#         # still counted as in hospital on day of discharge
#         start_1_date <= !!sym(v_discharged)   
#     )
#   
# }
# 
# data_start_1_date <- data_start_1_date %>%
#   transmute(
#     patient_id, day = 0, 
#     inhosp_planned, inhosp_unplanned, inhosp_covidunplanned,
#     # eol care initiated
#     endoflife = !is.na(endoflife_date) &
#       (endoflife_date < start_1_date),
#     # positive test in past 30 days
#     recentpostest = !is.na(postest_1_date) &
#       (postest_1_date < start_1_date) & 
#       as.integer(start_1_date - postest_1_date) < 30
#   )
# 
# ## each time status updates after start_1_date
# 
# static_vars <- c("patient_id", "start_1_date")
# 
# # end of life:
# data_endoflife <- extract %>%
#   select(all_of(static_vars), endoflife_date) %>%
#   # remove missing rows as already coded as 0 on start date
#   filter(!is.na(endoflife_date)) %>%
#   # remove endoflife_date < start_1_date as already coded as 1 on start date
#   filter(start_1_date <= endoflife_date) %>%
#   transmute(
#     patient_id, 
#     day = as.integer(endoflife_date - start_1_date),
#     endoflife = TRUE
#     )
# 
# # postest:
# # transform data to long format
# data_postest_long <- extract %>%
#   select(all_of(static_vars), starts_with("postest")) %>%
#   pivot_longer(cols = -all_of(static_vars), values_drop_na = TRUE) 
# 
# data_postest <- data_postest_long %>%
#   select(all_of(static_vars), postest_start_date = value) %>%
#   # end date for a "positive test episode" is 30 days after date of positive test
#   mutate(postest_end_date = postest_start_date+30) %>%
#   # get the date of the next positive test
#   group_by(patient_id) %>%
#   mutate(lead_postest_start_date = lead(postest_start_date)) %>%
#   ungroup() %>%
#   # if the next positive test is less than 30 days after the current one, remove the end date
#   mutate(
#     across(
#       postest_end_date, 
#       ~if_else(
#         !is.na(lead_postest_start_date) & (lead_postest_start_date <= .x),
#         as.Date(NA_character_),
#         .x
#         ))) %>%
#   select(-lead_postest_start_date) %>%
#   # transform data to long format
#   pivot_longer(cols = -all_of(static_vars), values_drop_na = TRUE) %>%
#   # postest=1 for start dates, =0 for end dates
#   mutate(postest = name == "postest_start_date") %>%
#   # arrange events (start and end dates) by their date, within individuals
#   group_by(patient_id) %>%
#   arrange(value, .by_group = TRUE) %>%
#   # calculate a "moving sum" of the current row and the previous row, 
#   # to flag instances where there are two start dates in a row
#   mutate(cumulsum = stats::filter(postest, rep(1,2), sides = 1)) %>%
#   ungroup() %>% 
#   # remove start dates if preceding row is also a start date (within patients)
#   filter(is.na(cumulsum) | cumulsum == 1) %>%
#   transmute(
#     patient_id,
#     day = as.integer(value - start_1_date),
#     postest
#   )
# 
# # check that there's never a new positive test less than 30 days after the last one
# cat("Check all postitive tests separated by at least 30 days.\n")
# cat("Following dataset should have 0 rows:\n")
# data_postest %>%
#   filter(postest) %>%
#   group_by(patient_id) %>%
#   mutate(timesincelast_postest = c(NA_integer_, diff(day))) %>%
#   ungroup() %>%
#   filter(timesincelast_postest < 30) %>% 
#   nrow()
# 
# # hospitalisations
# data_hosp <- extract %>%
#   select(all_of(static_vars), matches("\\w+_planned_\\d+_date")) %>%
#   pivot_longer(
#     cols = -all_of(static_vars),
#     values_drop_na = TRUE
#   ) %>%
#   # remove rows before start date as already coded
#   filter(start_1_date <= value) %>%
#   mutate(day = as.integer(value - start_1_date)) %>%
#   transmute(
#     patient_id, 
#     # still counted as in hospital on day of discharge
#     day = if_else(str_detect(name, "^admitted"), day, day+1),
#     inhosp_planned = str_detect(name, "^admitted")
#   )
# 
# 
# data_fup <- data_endoflife %>%
#   full_join(data_postest, by = c("patient_id", "day")) %>%
#   full_join(data_hosp, by = c("patient_id", "day")) 
# 
# 
# bind_rows(
#   data_start_1_date,
#   data_fup
# ) %>% filter(patient_id<10) %>% arrange(patient_id, day) %>%
#   mutate(across(-c(patient_id, day), ~if_else(is.na(.x) & day == 0, FALSE, .x))) %>%
#   group_by(patient_id) %>%
#   fill(everything(), .direction = "down")

########
  
static_vars <- c("patient_id", "start_1_date")

# endoflife
# end of life:
data_endoflife <- extract %>%
  select(all_of(static_vars), endoflife_date) %>%
  transmute(
    patient_id, 
    endoflife = !is.na(endoflife_date),
    day = if_else(
      endoflife,
      as.integer(endoflife_date - start_1_date),
      0
    )
  )

# postest:
data_postest <- extract %>%
  select(all_of(static_vars), starts_with("postest")) %>%
  pivot_longer(cols = -all_of(static_vars), values_drop_na = TRUE) %>%
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
  # remove start dates if preceding row is also a start date (within patients)
  filter(is.na(cumulsum) | cumulsum == 1) %>%
  transmute(
    patient_id,
    day = as.integer(value - start_1_date),
    postest
  )
  
## clean hospitalisation data:
# issue1: admission_*_x_date > discharge_*_x_date: 
#         stay of length 1 starting on admission date
# issue2: admission_*_x_date == discharge_*_x_date: 
#         stay of length 1
# issue3: admission_*_x_date <= discharge_*_{x-1}_date: 
#         stay continuous lasting from admission_*_x_date until discharge_*_{x+1}_date

hosp_vars <- recurring_vars[str_detect(recurring_vars, "admitted")]
issues <- list()
data_hosp_processed <- extract %>% select(all_of(static_vars), starts_with(c("admitted", "discharged")))

for (v in hosp_vars) {
  
  indices <- data_hosp_processed %>% select(starts_with(v)) %>% names() %>% 
    str_extract(., "\\d+") %>% as.integer() %>% sort()
  
  v_type <- str_remove(v, "admitted_")
  v_discharged <- str_replace(v, "admitted", "discharged")
  
  issues[[v]] <- list(issue1 = integer(), issue2 = integer(), issue3 = integer())
  
  # summarise issues
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
  # way without excluding patients with this issue
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


# derive data_hosp
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
      !!sym(gluec("inhosp_{v_type}")) := str_detect(name, "^admitted")
    )
  
}

by_vars <- c("patient_id", "day")

data_all <- data_endoflife %>%
  full_join(data_postest, by = by_vars) %>%
  full_join(data_hosp$admitted_planned, by = by_vars) %>%
  full_join(data_hosp$admitted_unplanned, by = by_vars) %>%
  full_join(data_hosp$admitted_covidunplanned, by = by_vars)

cat("check for days with >1 rows:\n")
data_all %>% 
  group_by(patient_id, day) %>%
  mutate(n=n()) %>% 
  filter(n>1) %>%
  nrow()

data_filled <- data_all %>%
  select(patient_id, day, everything()) %>%
  arrange(patient_id, day) %>%
  group_by(patient_id) %>%
  mutate(sequence = row_number()) %>%
  ungroup() %>%
  mutate(
    across(
      -c(patient_id, day), 
      ~if_else((sequence == 1) & (is.na(.x)), FALSE, .x)
      )
    ) %>%
  select(-sequence) %>%
  group_by(patient_id) %>%
  fill(-c(patient_id, day), .direction = "down")

data_0 <- data_filled %>%
  filter(day <= 0) %>%
  group_by(patient_id) %>%
  summarise(across(everything(), last)) %>%
  mutate(day=0)

data_final <- bind_rows(
  data_0,
  data_filled %>% filter(day > 0)
) %>%
  arrange(patient_id, day) 


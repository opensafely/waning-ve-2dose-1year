library(tidyverse)
library(rlang)


# preliminaries ----
study_parameters <- readr::read_rds(here::here("analysis", "lib", "study_parameters.rds"))

set.seed(study_parameters$seed)

n <- study_parameters$n

K <- study_parameters$K

start_date <- study_parameters$start_date

regions <- readr::read_csv(here::here("analysis", "lib", "regions.csv"))

# function for binary flag variables
var_binary <- function(.data, name, incidence=0.05, keep_vars = TRUE) {
  out <- .data %>%
    mutate({{name}} := as.integer(rbernoulli(n=nrow(.), p=incidence)))
  
  if (! keep_vars) {
    out <- out %>% select({{name}})
  }
  
  return(out)
  
}

# function for date variables
var_date <- function(.data, name, incidence=0.05, 
                     earliest="2000-01-01", latest="2021-12-31",
                     keep_vars = TRUE) {
  
  # check date in correct format
  d <- lapply(list(earliest, latest), function(x) try(as.Date(x, format="%Y-%m-%d")))
  if(any("try-error" %in% sapply(d, class)) || any(sapply(d, is.na))) {
    stop("`earliest` and `latest` must be in %Y-%m-%d format.")
  } 
  
  out <- .data %>%
    mutate({{name}} := if_else(
      rbernoulli(n=nrow(.), p=incidence), 
      sample(
        x = seq(as.Date(earliest), as.Date(latest), 1), 
        size = nrow(.),
        replace = TRUE),
      NA_Date_))
  
  if (! keep_vars) {
    out <- out %>% select({{name}})
  }
  
  return(out)
}

# function for categorical variables
var_category <- function(.data, name, categories, conditions=NULL, ratios=NULL, keep_vars = TRUE) {

  if (is.null(conditions)) {

    if (is.null(ratios)) ratios <- rep(1/length(categories), length(categories))

    out <- .data %>%
      mutate({{name}} := sample(x = categories, size = nrow(.), replace=TRUE, prob = ratios))

  } else {

    patterns_string <- str_c(str_c(conditions, "~\"", categories, "\""), collapse = "; ")

    patterns <- parse_exprs(patterns_string)

    out <- .data %>%
      mutate({{name}} := case_when(!!!patterns))

  }
  
  if (! keep_vars) {
    out <- out %>% select({{name}})
  }
  
  return(out)

}

# function for recurrent bmi variables (where 1 is most recent, 2 2nd most recent etc.)
vars_bmi_recurrent <- function(.data, r = 6, decay = 0.2) {
  
  out <- .data %>%
    var_date(bmi_0_date_measured, 
             incidence = 0.95, 
             earliest = "2018-01-01",
             keep_vars = FALSE) %>%
    mutate(
      bmi_0 = if_else(
        is.na(bmi_0_date_measured),NA_real_, rnorm(nrow(.), 28, 8)),
      bmi_0_stage = sample(
        x = c(NA_character_, str_c("Obese class ", c("I", "II", "III"))),
        size = nrow(.),
        replace=TRUE,
        prob = c(0.7,0.1,0.1,0.1)),
      bmi_0_stage_date = if_else(
        is.na(bmi_0_stage),
        as.Date(NA_character_),
        as.Date("2021-01-01")
      )
    )
  
  for (i in 1:r) {
    out <- out %>%
      mutate(missing = sample(
        x=c(NA_real_, 1),
        size = nrow(.), 
        replace=TRUE, 
        prob = c(decay, 1-decay)
        )) %>%
      mutate(
        !! glue("bmi_{i}_date_measured") := 
          !! sym(glue("bmi_{i-1}_date_measured")) 
        + sample(x=1:500, size = nrow(.), replace = TRUE)*missing,
        !! glue("bmi_{i}") := !! sym(glue("bmi_{i-1}"))
        + rnorm(n = nrow(.))*missing,
        !! glue("bmi_{i}_stage") := !! sym(glue("bmi_{i-1}_stage")),
        !! glue("bmi_{i}_stage_date") := !! sym(glue("bmi_{i-1}_stage_date"))
        ) 
  }
  
  return(out %>% select(-missing))
  
}

# recurrent date variable, where var_1 is the earliest date
var_date_recurrent <- function(
  .data, 
  name_string,
  incidence, 
  earliest = start_date, 
  r = 5, 
  decay = 0.5,
  start_index = 0) {
  
  out <- .data %>%
    var_date(
      name = !! sym(glue("{name_string}_{start_index}_date")),
      incidence = incidence, 
      earliest = earliest,
      keep_vars = FALSE)
  
  for (i in (start_index+1):r) {
    out <- out %>%
      mutate(missing = sample(
        x=c(NA_real_, 1),
        size = nrow(.), 
        replace=TRUE, 
        prob = c(decay, 1-decay)
      )) %>%
      mutate(
        !! glue("{name_string}_{i}_date") := !! sym(glue("{name_string}_{i-1}_date"))
        + sample(x=1:50, size = nrow(.), replace = TRUE)*missing
      ) 
  }
  
  return(out %>% select(-missing))
}

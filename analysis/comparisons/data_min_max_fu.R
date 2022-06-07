################################################################################

# This script:
# - calculates min and max follow-up dates for each subgroup

################################################################################

library(tidyverse)

################################################################################
# create output directories ----
fs::dir_create(here::here("output", "lib"))

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

################################################################################
## read data
# covariates data
data_all <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) 

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
data_min_max_fu <- data_all %>%
  group_by(subgroup) %>%
  summarise(
    min_fu_date = min(start_1_date),
    max_fu_date = max(end_6_date),
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

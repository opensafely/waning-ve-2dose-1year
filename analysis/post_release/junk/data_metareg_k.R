library(tidyverse)

################################################################################
if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

################################################################################

metareg_results_0 <- haven::read_dta(here::here(release_folder, "results.dta"))

metareg_results <- metareg_results_0 %>%
  rename(subgroup = stratum, comparison = vaccine) %>%
  mutate(across(subgroup,
                factor,
                levels = 1:4,
                labels = subgroups)) %>%
  mutate(across(sex,
                factor,
                levels = 1:3,
                labels = c("Both", "Female", "Male"))) %>%
  mutate(across(ageband,
                factor,
                levels = 1:3,
                labels = c("all", "65-74 years", "75+ years"))) %>%
  mutate(across(comparison, 
                factor,
                levels = 1:3,
                labels = c("BNT162b2", "ChAdOx1", "both")
                )) %>%
  mutate(across(outcome,
                factor,
                levels = 1:5,
                labels = c("covidadmitted",
                           "coviddeath",
                           "postest",
                           "noncoviddeath",
                           "anytest"))) %>%
  mutate(
    # slope ci
    logrhr_lower = logrhr - qnorm(0.975)*selogrhr,
    logrhr_higher = logrhr + qnorm(0.975)*selogrhr,
    # intercept ci
    loghr1_lower = loghr1 - qnorm(0.975)*seloghr1,
    loghr1_higher = loghr1 + qnorm(0.975)*seloghr1,
    ) 

metareg_results_k <- metareg_results %>%
  distinct(subgroup, sex, ageband, comparison, outcome) %>%
  mutate(k=factor(1, levels=1:6)) %>%
  complete(subgroup, sex, ageband, comparison, outcome, k) %>%
  left_join(
    metareg_results,
    by = c("subgroup", "sex", "ageband", "comparison", "outcome")
  ) %>%
  mutate(across(k, ~as.integer(as.character(.x)))) %>%
  mutate(
    line = loghr1 + (k-1)*logrhr#, # use k-1 because intercept at k=-1
    # CI for slope not valid as don't know covariance between intercept and slope
    # line_lower = loghr1_lower + (k-1)*logrhr_lower,
    # line_higher = loghr1_higher + (k-1)*logrhr_higher
  ) %>%
  mutate(across(starts_with("line"), exp))

readr::write_rds(
  metareg_results_k,
  here::here(release_folder, "metareg_results_k.rds")
)

metareg_results_rhr <- metareg_results %>%
  distinct(outcome, subgroup, sex, ageband, comparison, logrhr, logrhr_lower, logrhr_higher) %>%
  mutate(across(starts_with("log"), exp)) %>%
  rename_with(~str_remove(.x, "log"))

readr::write_rds(
  metareg_results_rhr,
  here::here(release_folder, "metareg_results_rhr.rds")
)


library(tidyverse)
library(glue)

# define release folder
release_folder <- here::here("release20220622")

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

###############################################################################
# cumulative incidence data
survtable_redacted <- readr::read_csv(
  here::here(release_folder, "survtable_redacted.csv")) %>%
  mutate(across(time, ~.x+2)) %>%
  mutate(across(c.inc, ~100*round(.x,2))) %>%
  mutate(across(subgroup,
                ~ case_when(
                  str_detect(.x, "65") ~ 1,
                  str_detect(.x, "18-64") ~ 2,
                  str_detect(.x, "40-64") ~ 3,
                  str_detect(.x, "18-39") ~ 4,
                  TRUE ~ NA_real_
                ))) 

survtable_redacted %>%
  filter(
    subgroup == 2,
    # arm == "BNT162b2",
    arm == "ChAdOx1",
    # arm == "Unvaccinated",
    # time>45
    time>20, time <40
  ) %>% print(n=50)

###############################################################################
# event counts
event_counts <- readr::read_csv(here::here(release_folder, "event_counts_all.csv")) 

event_counts %>%
  filter(
    subgroup == 2,
    outcome == "coviddeath"
  ) %>%
  select(-outcome) %>%
  arrange(arm, k) %>%
  print(n=Inf)



# read estimates data
estimates_all <- readr::read_csv(here::here(release_folder, "estimates_all.csv")) 


estimates_all %>%
  filter(
    subgroup == 1,
    comparison != "both",
    model == "adjusted",
    variable == "k",
    !reference_row
    ) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~format(round(exp(.x), 2), nsmall=2))) %>%
  arrange(outcome, comparison, period) %>%
  transmute(
    comparison, outcome, period, n_obs_model,
    value = glue("{estimate} ({conf.low}, {conf.high})")
  ) %>%
  print(n=Inf)

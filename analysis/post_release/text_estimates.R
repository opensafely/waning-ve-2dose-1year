################################################################################
## setup 
library(tidyverse)

release_folder <- here::here("release20221006")

################################################################################
## read and derive metadata 
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

# labels for comparison periods
K <- study_parameters$K
ends <- seq(2, (K+1)*4, 4)
starts <- ends + 1
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")

################################################################################
## read and process results 
# read event counts
event_counts <- readr::read_csv(
  here::here(release_folder, "event_counts_all.csv")
)

# read subsequent vax
survtable_redacted <- readr::read_csv(
  here::here(release_folder, "survtable_redacted.csv")
) %>%
  # the following line to rescale to weeks after second dose rather than weeks after start_1_date
  mutate(across(time, ~.x+2))

# read all hr estimates and filter estimates of vaccine effectiveness
estimates_all <- readr::read_csv(
  here::here(release_folder, "estimates_all.csv")
) %>%
  mutate(across(subgroup, factor, 
                levels = seq_along(subgroups),
                labels = subgroups)) %>%
  filter(
    variable == "k", # comparison period estimates
    !reference_row, 
    model=="adjusted", # keep only max adjusted model
    outcome %in% c("postest", "covidadmitted", "coviddeath", "noncoviddeath")
  ) %>%
  # calculate exponential of estimate and confidence intervals
  mutate(across(c(estimate, conf.low, conf.high), ~round(exp(.x),3))) %>%
  mutate(across(label, ~factor(as.numeric(.x), labels = weeks_since_2nd_vax)))

################################################################################
## study population subsection 
# number of people in each vaccine group
event_counts %>%
  filter(k==1, outcome=="coviddeath") %>%
  group_by(arm) %>%
  summarise(sum(n)) %>%
  ungroup()

# maximum follow-up in each subgroup following second dose
survtable_redacted %>%
  filter(arm!="Unvaccinated") %>%
  group_by(subgroup) %>%
  summarise(max(time))

################################################################################
## cumulative incidence of subsequent dose subsection
# maximum number of weeks each subgroup is followed up for following 2nd dose
survtable_redacted %>%
  filter(arm=="Unvaccinated") %>%
  group_by(subgroup) %>%
  arrange(desc(time), .by_group=TRUE) %>% 
  slice(1:5) %>%
  mutate(across(c.inc, ~.x*100))

# summary of third dose update (alter subgroup and time limits accordingly)
survtable_redacted %>%
  filter(
    # str_detect(subgroup,"18-39"), 
    arm!="Unvaccinated",
    time > 26 & time < 35
    ) %>%
  mutate(across(c.inc, ~.x*100)) %>%
  select(subgroup, arm, time, c.inc) %>%
  pivot_wider(names_from=arm, values_from=c.inc) %>%
  print(n=100)
  
################################################################################
## waning vaccine effectiveness

estimates_all %>%
  filter(
    as.integer(subgroup)==4,
    outcome%in%c("covidadmitted", "coviddeath", "postest", "noncoviddeath")
    ) %>%
  transmute(subgroup=as.integer(subgroup), comparison, outcome, label, 
            value = glue::glue("{format(round(estimate, 2), nsmall=2)} ({format(round(conf.low, 2), nsmall=2)}, {format(round(conf.high, 2), nsmall=2)})")) %>%
  pivot_wider(names_from=comparison, values_from = value) %>%
  arrange(outcome,label) %>%
  print(n=Inf)

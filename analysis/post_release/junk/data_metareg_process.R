library(tidyverse)

if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")

# read data
estimates_all <- bind_rows(
  readr::read_csv(file.path(release_folder, "estimates_all.csv")),
  readr::read_csv(file.path(release_folder, "estimates_6575.csv"))
)

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds")
)
outcomes <- outcomes[outcomes!="covidemergency"]

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

# process
data_metareg_0 <- estimates_all %>%
  filter(
    outcome %in% outcomes,
    !reference_row,
    variable %in% "k",
    model %in% "adjusted" # error in the labeling, this corresponds to unadjusted
    ) %>%
  mutate(
    model = "adjusted",
    sex = if_else(
      str_detect(subgroup, "Female|Male"),
      str_extract(subgroup, "Female|Male"),
      "Both"
      ),
    ageband = case_when(
      str_detect(subgroup, "_65") ~ "65-74 years",
      str_detect(subgroup, "_75") ~ "75+ years",
      TRUE ~ "all"
    ),
    subgroup = as.integer(str_extract(subgroup, "\\d"))
    ) %>%
  select(subgroup, sex, ageband, comparison, outcome, model, k = label, estimate, conf.low, conf.high) %>%
  mutate(across(subgroup, factor, levels = 1:4, labels = subgroups)) 

# expand to all combinations
data_metareg <- data_metareg_0 %>%
  expand(subgroup, comparison, sex, ageband, outcome, model, k) %>%
  filter(!(subgroup %in% "18-39 years" & comparison != "BNT162b2")) %>%
  left_join(data_metareg_0) %>%
  # following lines because in a few instances seloghr=0, and this is breaking the metareg code 
  mutate(seloghr = (conf.high - conf.low)/(2*1.96)) %>%
  mutate(across(c(estimate, conf.low, conf.high),
                ~if_else(near(seloghr, 0), NA_real_, .x))) %>%
  select(-seloghr)

# save data
readr::write_csv(
  data_metareg,
  here::here(release_folder, "data_metareg.csv"))

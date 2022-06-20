library(tidyverse)
library(glue)
library(lubridate)
library(RColorBrewer)

cat("setup\n")

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  subgroup_label <- 1L
  
} else{
  subgroup_label <- as.integer(args[[1]])
}

################################################################################
# output directories
fs::dir_create(here::here("output", "tte", "images"))

################################################################################
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

# function to be applied in dplyr::filter
occurs_after_start_date <- function(cov_date, index_date) {
  is.na(cov_date) | index_date < cov_date
}

################################################################################
cat("read and process data\n")

data_tte <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) %>%
  select(patient_id, subgroup, arm, split, death_date, dereg_date, subsequent_vax_date,
         starts_with(c("start", "end"))) %>%
  filter(subgroup %in% subgroups[subgroup_label]) %>%
  pivot_longer(
    cols = c(starts_with(c("start","end"))),
    names_to = c(".value", "k"),
    names_pattern = "(.*)_(.*)_date"
    ) %>%
  mutate(across(k, as.integer)) %>%
  filter(
    is.na(split) |
      ((k %% 2) == 0 & split == "even") |
      ((k %% 2) != 0 & split == "odd")
  ) %>%
  select(-split) %>%
  filter_at(
    vars(str_c(unique(c("subsequent_vax", "dereg", "death")), "_date")),
    all_vars(occurs_after_start_date(cov_date = ., index_date = start))
  ) %>%
  filter(start <=  as.Date(study_parameters$end_date)) %>%
  mutate(across(ends_with("date"), ~if_else(start <= .x & .x <= end, .x, as.Date(NA_character_)))) %>%
  mutate(across(end, ~pmin(end, death_date, dereg_date, subsequent_vax_date, as.Date(study_parameters$end_date), na.rm = TRUE))) %>%
  select(patient_id, subgroup, arm, k, start, end) %>%
  mutate(across(k, factor, levels = 1:K))

################################################################################
cat("derive plot data\n")

cat("define variant dates\n")
delta_start <- as.Date("2021-06-01")
omicron_start <- as.Date("2021-12-01")
delta_end <- as.Date("2021-12-15")

cat("earliest and latest fu dates\n")
min_date <- min(data_tte$start)
max_date <- max(data_tte$start)

cat("prepare data\n")
cat("wide data with each individual's start and end date\n")
data_tte_wide <- data_tte %>%
  mutate(
    tstart = as.integer(start - min_date),
    tstop =  as.integer(end - min_date) 
  ) %>%
  select(-start, -end) 

cat("split into BNT162b2, ChAdOx1 and unvax\n")
# because it was failing on full dataset due to size when one row per person-day
data_patients_list <- data_tte_wide %>% 
  group_split(arm)

data_tte_long_list <- list()
for (i in seq_along(data_patients_list)) {
  
  cat(glue("{i}:\n"))
  cat("every possible follow-up day for each individual\n")
  data_patients <- data_patients_list[[i]] %>%
    distinct(patient_id) %>%
    uncount(weights = max(data_tte_wide$tstop)+1) %>%
    mutate(
      time = rep(
        0:max(data_tte_wide$tstop), 
        times = n_distinct(data_patients_list[[i]]$patient_id)
        )
      )
  
  cat("keep only follow-up days that are within start and end date for each individual\n")
  data_tte_long_list[[i]] <- data_tte_wide %>%
    left_join(
      data_patients, by = "patient_id"
    ) %>%
    filter(tstart < time, time <= tstop) %>%
    group_by(k, time) %>%
    count() %>%
    ungroup() %>%
    mutate(date = min_date + days(time)) %>%
    select(-time) %>%
    # round up to nearest 7
    # then divide by 1000 to keep the plot tidy (will add note to y-axis)
    mutate(across(n, ~ceiling_any(.x, to = 7)/1000)) 
  
}

cat("bind across i\n")
data_tte_long <- bind_rows(data_tte_long_list)

# save data
write_csv(
  data_tte_long,
  here::here("output", "tte", "images", glue("check_fu_{subgroup_label}.csv"))
)

################################################################################
cat("metadata and objects for plot\n")
col_palette <- brewer.pal(n = 3, name = "Dark2")
xintercepts <- Date()
names_xintercepts <- character()
n_mult <- numeric()
k_print <- integer()
index <- integer()
if (min_date <=  delta_start) {
  xintercepts <- c(xintercepts, delta_start)
  names_xintercepts <- "Delta became\ndominant"
  n_mult <- c(n_mult, 0.75)
  k_print <- c(k_print, 12)
  index <- c(index, 1)
}
if (omicron_start <=  max_date) {
  xintercepts <- c(xintercepts, omicron_start)
  names_xintercepts <- c(names_xintercepts, "First cases of\nOmicron detected")
  n_mult <- c(n_mult, 0.75)
  k_print <- c(k_print, 11)
  index <- c(index, 2)
}
if (delta_end <=  max_date) {
  xintercepts <- c(xintercepts, delta_end)
  names_xintercepts <- c(names_xintercepts, "50% of cases\nOmicron")
  n_mult <- c(n_mult, 0.75)
  k_print <- c(k_print, 12)
  index <- c(index, 3)
}
names(xintercepts) <- names_xintercepts

ann_text <- tibble(
  date = xintercepts,
  n = n_mult*max(data_tte_long$n),
  k = factor(k_print, levels = 1:K),
  lab = names(xintercepts),
  lab_col = col_palette[1:length(xintercepts)]
)

cat("create plot\n")
p <- data_tte_long %>%
  ggplot(aes(x = date, y = n)) +
  geom_bar(stat = "identity", alpha = 0.5, width=1) +
  geom_vline(
    data = bind_rows(
      lapply(
        1:max(as.integer(as.character(data_tte_long$k))),
        function(x) 
          ann_text %>% mutate(k=x)
      )
    ) %>%
      mutate(across(k, factor)), 
    aes(xintercept = date, colour = lab_col),
    linetype = "dashed") +
  labs(
    x = "Date of follow-up", 
    y = "Number of individuals followed up (x 1000)"
  ) +
  facet_grid(k~.) +
  geom_label(
    data = ann_text, 
    aes(label = lab, colour = lab_col),
    size = 2
  ) +
  scale_color_manual(
    values = col_palette[index]
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%d %b %y",
    limits = c(min_date - days(14), max_date + days(14))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_text(
      size=10, 
      margin = margin(t = 10, r = 0, b = 0, l = 0)
    ),
    axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position = "none"
  )

ggsave(p,
       filename = here::here("output", "tte", "images", glue("check_fu_{subgroup_label}.png")),
       width=15, height=20, units="cm")
  
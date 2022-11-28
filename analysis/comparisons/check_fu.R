library(tidyverse)
library(glue)
library(lubridate)
library(RColorBrewer)

cat("setup\n")

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  release_folder <- here::here("release20221006")
  subgroup_label <- 1L
  
} else{
  subgroup_label <- as.integer(args[[1]])
}

################################################################################
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

cat("define variant dates\n")
delta_start <- as.Date("2021-06-01")
omicron_start <- as.Date("2021-12-01")
delta_end <- as.Date("2021-12-15")

# if running locally read extracted data:
if(Sys.getenv("OPENSAFELY_BACKEND") %in% "") {
  
  if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")
  image_path <- file.path(release_folder, "images")
  
  data_tte_long <- readr::read_csv(file.path(release_folder, glue("check_fu_{subgroup_label}.csv"))) %>%
    mutate(across(k, factor, levels = 1:K))
  
  
} else {
  # otherwise derive the data
  
  ################################################################################
  # output directories
  image_path <- here::here("output", "tte", "images")
  fs::dir_create(image_path)
  
  ################################################################################
  # redaction functions
  source(here::here("analysis", "functions", "redaction_functions.R"))
  
  # function to be applied in dplyr::filter
  occurs_after_start_date <- function(cov_date, index_date) {
    is.na(cov_date) | index_date < cov_date
  }
  
  ################################################################################
  cat("read and process data\n")
  
  # read start dates
  data_all <- readr::read_rds(
    here::here("output", "data", "data_all.rds")) %>%
    select(patient_id, start_1_date)
  
  # read tte data for BNT162b2 and unvax
  data_tte <-  readr::read_rds(
    here::here("output", "tte", "data", glue("data_tte_BNT162b2_{subgroup_label}_coviddeath.rds"))
  )
  
  # if subgroups 1:3, read tte data for ChAdOx1
  if (subgroup_label %in% 1:3) {
    data_tte <- bind_rows(
      data_tte,
      readr::read_rds(
        here::here("output", "tte", "data", glue("data_tte_ChAdOx1_{subgroup_label}_coviddeath.rds"))
      ) %>%
        filter(arm %in% "ChAdOx1")
    )
  }
  
  # going tte and start dates
  data_tte <- data_tte %>%
    left_join(data_all, by = "patient_id") %>%
    mutate(across(k, factor, levels = 1:K))
  
  ################################################################################
  cat("derive plot data\n")
  
  cat("earliest and latest fu dates\n")
  # calculate min_tstart just in case the person with tstart = 0 has been removed
  min_tstart <- min(data_tte$tstart) 
  # origin_date was used in data_tte to calculate tstart tstop so now use it back calculate dates
  origin_date <- min(data_tte$start_1_date) - days(min_tstart) 
  
  cat("split into BNT162b2, ChAdOx1 and unvax\n")
  # because it was failing on full dataset due to size when one row per person-day
  data_patients_list <- data_tte %>% 
    mutate(tmp_group = str_c(arm, k, sep="_")) %>%
    group_split(tmp_group)
  
  data_tte_long_list <- list()
  for (i in seq_along(data_patients_list)) {
    
    cat(glue("{i}:\n"))
    
    # prepare data
    data_tte_long_list[[i]] <- data_patients_list[[i]] %>%
      group_by(patient_id, k) %>%
      nest() %>%
      mutate(
        date = map(data, ~.x$tstart:.x$tstop)
      ) %>%
      unnest(cols = c(data, date)) %>%
      ungroup() %>%
      mutate(date = origin_date + days(date)) %>%
      group_by(k, date) %>%
      count() %>%
      ungroup() %>%
      # round up to nearest 7
      # also divide by 1000 to keep the plot tidy (will add note to y-axis)
      mutate(across(n, ~ceiling_any(.x, to = 7)/1000)) 
    
  }
  
  cat("bind across i\n")
  data_tte_long <- bind_rows(data_tte_long_list)
  
  # save data
  write_csv(
    data_tte_long,
    here::here("output", "tte", "images", glue("check_fu_{subgroup_label}.csv"))
  )
  
}

################################################################################
cat("metadata and objects for plot\n")
col_palette <- brewer.pal(n = 3, name = "Dark2")
xintercepts <- Date()
names_xintercepts <- character()
n_mult <- numeric()
k_print <- integer()
index <- integer()

cat("earliest and latest fu dates\n")
min_date <- min(data_tte_long$date)
max_date <- max(data_tte_long$date)


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
    limits = c(min_date, max_date)
    # limits = c(min_date - days(14), max_date + days(14))
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
       filename = file.path(image_path, glue("check_fu_{subgroup_label}.png")),
       width=15, height=20, units="cm")
  
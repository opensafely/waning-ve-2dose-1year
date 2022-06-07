################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)

################################################################################
fs::dir_create(here::here("output", "models_cox", "images"))

################################################################################
# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

################################################################################
estimates_all <- readr::read_csv(here::here("output", "release_objects", "estimates_all.csv")) %>%
  filter(variable != "k", !reference_row) %>%
  mutate(across(c(estimate, conf.low, conf.high), exp)) %>%
  mutate(
    variable = glue("{variable}_{label}"),
    sex = case_when(
      str_detect(subgroup, "Female") ~ "Female",
      str_detect(subgroup, "Male") ~ "Male",
      TRUE ~ "Both"
    ),
    subgroup = as.integer(str_extract(subgroup, "\\d"))
  ) %>%
  mutate(across(subgroup, factor, levels = subgroup_labels, labels = subgroups)) 


################################################################################
plot_coefs <- function(c, i, s) {
  
  plot_data <- estimates_all %>%
    filter(
      comparison %in% c,
      subgroup %in% subgroups[i],
      sex %in% s
    ) %>%
    mutate(
      comparison = glue("{comparison}: sex={s}, {subgroups[i]}"),
      short_term = variable,
      period = as.factor(period)
    ) %>%
    mutate(across(short_term, ~str_remove(.x, " \\w+ deprived"))) %>%
    mutate(across(short_term, ~str_remove(.x, "ethnicity"))) %>%
    mutate(across(short_term, ~str_remove(.x, "bmi"))) %>%
    mutate(across(short_term, ~str_remove(.x, "TRUE"))) %>%
    mutate(across(short_term, ~str_trunc(.x, width = 15, side = "center"))) 
  
  # plot
  plot_data %>%
    ggplot(aes(x = short_term, y = estimate, colour = period)) +
    geom_hline(yintercept = 1, colour = "black") +
    geom_point(alpha = 0.3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 1, alpha = 0.5) +
    facet_grid(outcome ~ comparison, scales = "free_x") +
    scale_y_log10(
      name = "hazard ratio",
      breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
      labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10"),
      # limits = c(min(0.1, y_min), max(10, y_max)),
      oob = scales::oob_keep
    ) +
    labs(
      x = NULL
      ) +
    guides(colour = guide_legend(nrow = 1)) +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, size = 8),
      axis.text.y = element_text(size = 6),
      legend.position = "bottom"
    )
  ggsave(filename = here::here("output", "models_cox", "images", glue("coefs_{c}_{i}_{s}.png")),
         width=20, height=30, units="cm")
  
}

for (c in c("BNT162b2", "ChAdOx1", "both")) {
  for (i in 1:4) {
    for (s in c("Both", "Male", "Female")) {
      try(plot_coefs(c=c,i=i,s=s))
    }
  }
}

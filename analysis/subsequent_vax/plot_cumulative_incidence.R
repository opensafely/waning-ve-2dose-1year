library(tidyverse)
library(glue)
library(survival)
library(lubridate)

################################################################################
## import study_parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)
subgroups_long <- if_else(
  subgroups %in% subgroups[c(3,4)],
  as.character(glue("{subgroups}*")),
  subgroups
)
subgroups_long_wrap <- str_wrap(subgroups_long, width = 25)

################################################################################

# if running locally read extracted data:
if(Sys.getenv("OPENSAFELY_BACKEND") %in% "") {

  if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")
  image_path <- file.path(release_folder, "images")

  survtable_redacted <- readr::read_csv(
    here::here(release_folder, "survtable_redacted.csv")) %>%
    mutate(across(subgroup,
                  ~ case_when(
                    str_detect(.x, "65") ~ 1,
                    str_detect(.x, "18-64") ~ 2,
                    str_detect(.x, "40-64") ~ 3,
                    str_detect(.x, "18-39") ~ 4,
                    TRUE ~ NA_real_
                  ))) %>%
    mutate(across(subgroup,
                  factor,
                  levels = subgroup_labels,
                  labels = subgroups_long_wrap))

} else { # else derive the data
  
  # create output directories
  fs::dir_create(here::here("output", "subsequent_vax", "images"))
  fs::dir_create(here::here("output", "subsequent_vax", "tables"))
  
  # redaction function for KM curves
  source(here::here("analysis", "functions", "round_km.R"))
  
  # read data
  data_all <- readr::read_rds(
    here::here("output", "data", "data_all.rds")) %>%
    rename("end_K_date" = glue("end_{K}_date")) %>%
    select(patient_id, subgroup, arm, start_1_date, end_K_date, 
           subsequent_vax_date, dereg_date, death_date)
  
  image_path <- here::here("output", "subsequent_vax", "images")
  
  # function to be applied in dplyr::filter
  no_evidence_of <- function(cov_date, index_date) {
    is.na(cov_date) | index_date < cov_date
  }
  
  data_tte <- data_all %>%
    mutate(across(subgroup,
                  factor,
                  levels = subgroups,
                  labels = subgroups_long_wrap)) %>%
    # mutate(
    #   # start date of comparison 1 
    #   start_fu_date = start_1_date,
    #   # end date of final comparison or end of data availability
    #   # end_fu_date = pmin(end_K_date, study_parameters$end_date) # already done in data_covariates_process
    # ) %>%
    # select(-start_1_date, -end_K_date) %>%
    rename(
      start_fu_date = start_1_date,
      end_fu_date = end_K_date
    ) %>%
    # remove if subsequent vaccine, death or dereg on or before start_fu_date
    filter_at(
      all_of(c("subsequent_vax_date", "death_date", "dereg_date")),
      all_vars(no_evidence_of(., start_fu_date))) %>%
    group_by(start_fu_date) %>%
    mutate(across(c(end_fu_date, subsequent_vax_date, death_date, dereg_date),
                  # time in weeks between start_fu_date and event
                  ~ as.integer(.x - start_fu_date)/7)) %>%
    ungroup() %>%
    select(-start_fu_date) %>%
    rename_at(vars(ends_with("_date")),
              ~ str_remove(.x, "_date")) %>%
    mutate(
      tte = pmin(subsequent_vax, dereg, death, end_fu, na.rm = TRUE),
      status = if_else(
        # because subsequent_vax must occur before death and dereg if occuring on same day
        !is.na(subsequent_vax) & subsequent_vax == tte,
        TRUE,
        FALSE
      )) %>%
    select(patient_id, arm, subgroup, tte, status) %>%
    mutate(across(arm, 
                  ~if_else(
                    arm == "unvax",
                    "Unvaccinated",
                    .x)
    ))
  
  survtable_redacted <- round_km(
    data = data_tte,
    time = "tte",
    event = "status",
    strata = c("subgroup", "arm"),
    threshold = 7
  ) %>%
    ungroup() %>%
    mutate(
      c.inc = 1-surv,
      c.inc.ll = 1-surv.ul,
      c.inc.ul = 1-surv.ll
    ) %>%
    mutate(across(starts_with("c."), round, digits = 5)) %>%
    select(subgroup, arm, time, n.risk, n.event, n.censor, c.inc) 
  
  readr::write_csv(
    survtable_redacted,
    here::here("output", "subsequent_vax", "tables", "survtable_redacted.csv"))
  
}

################################################################################
# scale for x-axis
x_breaks <- seq(0,K*4,4)
x_labels <- x_breaks + 2

# create plot
plot_out <- survtable_redacted %>%
  ggplot(aes(x = time, y = c.inc, colour = subgroup)) +
  geom_step(size=0.8) +
  facet_wrap(~ arm, nrow=2) +
  scale_color_viridis_d(
    name = "Subgroup"
  ) +
  scale_x_continuous(
    breaks = x_breaks, # scale is time since start of period 1
    labels = x_labels # label scale as time since second vax
  ) +
  labs(
    x = "Weeks since second dose",
    y = "Incidence of 1st dose                 Incidence of 3rd dose"
    ) +
  theme_bw() +
  theme(
    
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=8),
    
    axis.title.x = element_text(size = 10, margin = margin(t = 20, r = 0, b = 10, l = 0)),
    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0, size = 11),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = c(0.75,0.2),
    legend.key.size = unit(0.8, "cm")
    
  )

ggsave(plot = plot_out,
       filename = file.path(image_path, "ci_vax.png"),
       width=16, height=12, units="cm")
  
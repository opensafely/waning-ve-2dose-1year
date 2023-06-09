library(tidyverse)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

tidy_dose3model <- read_csv(here::here("release20230605", "tidy_dose3model.csv")) 

historyof_vars <-c(
  "Learning disability" = "learndis", 
  "Severe mental illness" = "sev_mental",
  "CD respiratory" = "crd",
  "CD heart" = "chd",
  "CD liver" = "cld",
  "CD kidney" = "ckd", 
  "CD neurological" = "cns", 
  "Diabetes" = "diabetes", 
  "Immunosuppression" = "immunosuppressed", 
  "Cancer" = "status_cancer",
  "Flu vaccine" = "flu_vaccine"
)

baseline_vars <- c(
  "Age*" = "age",
  "Sex" = "sex",
  "IMD Quintile" = "imd", 
  "Ethnicity" = "ethnicity", 
  "BMI" = "bmi", 
  "Brand of second dose" = "vax2_brand",
  "Test history**" =  "test_hist_n"
  )           

timeupdating_vars <- c(
  "End of life care initiated" = "status_endoflife", 
  "Positive SARS-CoV-2 test\nin past 30 days" = "status_postest",
  "In hospital on planned\nadmission" = "status_planned", 
  "In hospital on unplanned\nadmission" = "status_unplanned", 
  "In hospital on unplanned\nadmission with COVID-19" = "status_covidunplanned" 
  )

plot_vars <- timeupdating_vars


plot_data <- tidy_dose3model %>%
  filter(!str_detect(variable, "^ns")) %>%
  mutate(
    subgroup = factor(
      model,
      levels = subgroup_labels,
      labels = str_wrap(subgroups, width = 25),
      )
    ) %>%
  mutate(across(variable, ~if_else(str_detect(.x, "age_smc"), "age", .x))) %>%
  mutate(across(label, ~if_else(variable==.x | str_detect(.x, "TRUE$"), "", .x))) %>%
  mutate(across(label, ~if_else(.x=="age_smc", "Degree=1", .x))) %>%
  mutate(across(label, ~if_else(.x=="age_smc_squared", "Degree=2", .x))) %>%
  filter(variable %in% plot_vars) %>%
  mutate(across(variable, ~factor(.x, levels = unname(all_vars), labels = names(all_vars)))) %>%
  mutate(
    variable_label = fct_rev(
      fct_inorder(
        str_replace(paste(variable, label, sep=": "),  "\\:\\s$", "")
      )
    )
  ) %>%
  droplevels()

position_dodge_val <- 1

if ("imd" %in% plot_vars) {
  plot_data <- plot_data %>%
    # reorder the factor levels so IMD printed 1:5
    mutate(
      variable_label = fct_relevel(
        plot_data$variable_label, 
        "IMD Quintile: 4", "IMD Quintile: 3", "IMD Quintile: 2", 
        after = which(levels(plot_data$variable_label) == "IMD Quintile: 5 least deprived")
      )
    )
}

plot_data %>%
  ggplot(aes(x = estimate, y = variable_label, colour = subgroup)) +
  geom_vline(xintercept = 1, colour='grey') +
  geom_linerange(
    aes(xmin = conf.low, xmax = conf.high),
    position = position_dodge(width = position_dodge_val)
  ) +
  geom_point(
    aes(), 
    position = position_dodge(width = position_dodge_val),
    alpha = 0.5,
    size = 1
  ) +
  facet_grid(rows = "variable", scales = "free_y", space = "free_y") +
  scale_color_discrete(name = NULL) +
  scale_x_log10(name = "Hazard Ratio", limits = c(0.1,3)) +
  scale_y_discrete(name=NULL) +
  guides(colour = guide_legend(nrow=1, byrow = TRUE)) +
  theme_bw() +
  theme(
    strip.text = element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.ticks = element_blank(),
    plot.margin = margin(b=50),
    legend.position = c(0.3,-0.1)
  )

ggsave(
  filename = here::here("release20230605", "testplot.png"),
  width = 14, height = 20, units = "cm"
)
  

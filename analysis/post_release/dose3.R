library(tidyverse)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

release_dir <- "release20230613"

### plots of categorical covariates
tidy_dose3model <- read_csv(here::here(release_dir, "tidy_dose3model.csv")) 

create_plot <- function(vars_group) {
  
  if (vars_group == "historyof") {
    plot_vars <-c(
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
  }
  if (vars_group == "baseline") {
    plot_vars <- c(
      "Sex" = "sex",
      "IMD Quintile" = "imd", 
      "Ethnicity" = "ethnicity", 
      "BMI" = "bmi", 
      "Brand of second dose" = "vax2_brand",
      "Test history**" =  "test_hist_n"
    )           
  }
  if (vars_group == "timeupdating") {
    plot_vars <- c(
      "End of life care initiated" = "status_endoflife", 
      "Positive SARS-CoV-2 test\nin past 30 days" = "status_postest",
      "In hospital on planned\nadmission" = "status_planned", 
      "In hospital on unplanned\nadmission" = "status_unplanned", 
      "In hospital on unplanned\nadmission with COVID-19" = "status_covidunplanned" 
    )  
  }
  
  plot_data <- tidy_dose3model %>%
    filter(!str_detect(variable, "^ns|^poly")) %>%
    mutate(
      subgroup = factor(
        model,
        levels = subgroup_labels,
        labels = str_wrap(subgroups, width = 25),
      )
    ) %>%
    mutate(across(label, ~if_else(variable==.x | str_detect(.x, "TRUE$"), "", .x))) %>%
    mutate(across(reference_row, ~replace_na(.x, FALSE))) %>%
    filter(variable %in% plot_vars) %>%
    mutate(across(variable, ~factor(.x, levels = unname(plot_vars), labels = names(plot_vars)))) %>%
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
  if ("bmi" %in% plot_vars) {
    plot_data <- plot_data %>%
      # reorder the factor levels so IMD printed 1:5
      mutate(
        variable_label = fct_relevel(
          plot_data$variable_label, 
          "BMI: Obese II (35-39.9)", "BMI: Obese I (30-34.9)", "BMI: Not obese", 
          after = which(levels(plot_data$variable_label) == "BMI: Obese III (40+)")
        )
      )
  }
  
  p <- plot_data %>%
    ggplot(aes(x = estimate, y = variable_label, colour = subgroup)) +
    geom_vline(xintercept = 1, colour='grey') +
    geom_linerange(
      aes(xmin = conf.low, xmax = conf.high),
      position = position_dodge(width = position_dodge_val)
    ) +
    geom_point(
      aes(shape = reference_row), 
      position = position_dodge(width = position_dodge_val),
      size = 1
    ) +
    facet_grid(rows = "variable", scales = "free_y", space = "free_y") +
    scale_color_discrete(name = NULL) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1), guide = "none") +
    scale_x_log10(name = "Hazard Ratio (HR > 1: higher rate of third dose)", limits = c(0.1,3)) +
    scale_y_discrete(name=NULL) +
    guides(
      colour = guide_legend(nrow=1, byrow = TRUE)
      ) +
    theme_bw() +
    theme(
      strip.text = element_blank(),
      panel.border = element_blank(),
      panel.spacing = unit(1, "lines"),
      
      axis.title.x = element_text(margin = margin(t=10)),
      
      axis.ticks = element_blank(),
      
      plot.margin = margin(b=50),
      legend.position = c(0.3,-0.12)
    )
  
  ggsave(
    plot = p,
    filename = here::here(release_dir, glue::glue("{vars_group}.png")),
    width = 14, height = 20, units = "cm"
  )
  
  return(p)
  
}


create_plot("baseline")
create_plot("historyof")
create_plot("timeupdating")


### plots of continuous covariates
hrs_age <- read_csv(here::here(release_dir, "hrs_age.csv")) 



plot_data <- hrs_age %>%
  # left_join(
  #   hrs_age %>% filter(y_trans==1) %>% select(subgroup, median_age = x),
  #   by = "subgroup"
  # ) %>%
  mutate(across(subgroup, ~factor(.x, levels = seq_along(subgroups), labels = subgroups))) %>%
  group_by(subgroup) %>%
  mutate(
    y_max = max(y),
    x_at_max_y = mean(if_else(y==y_max, x, NA_real_), na.rm=TRUE)
  ) %>%
  # mutate(min_age = min(x)) %>%
  # mutate(y_min_age = mean(if_else(x==min_age,y,NA_real_), na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(
    y_trans = y,
    conf.ll = y-(1.96*se),
    conf.ul = y+(1.96*se)
  ) %>%
  # rescale y relative to value at median_val
  mutate(across(c(y_trans, conf.ll, conf.ul), ~exp(.x - y_max))) 


plot_data %>%
  ggplot(aes(x = x)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_ribbon(aes(ymin = conf.ll, ymax = conf.ul, fill = subgroup), alpha = 0.2) +
  geom_line(aes(y = y_trans, colour = subgroup)) +
  # geom_vline(aes(xintercept = median_age, colour = subgroup), linetype = "dashed") +
  # geom_vline(aes(xintercept = min_age, colour = subgroup), linetype = "dashed") +
  geom_vline(aes(xintercept = x_at_max_y, colour = subgroup), linetype = "dashed") +
  geom_label(aes(x = x_at_max_y, y = 0.6, label = x_at_max_y, colour = subgroup)) +
  facet_wrap(facets = "subgroup", scales = "free_x") +
  scale_y_log10() +
  labs(x = "Age in years", y = "Hazard ratio") +
  scale_color_discrete(guide = "none") +
  scale_fill_discrete(guide = "none") +
  theme_bw() +
  theme(
    axis.title.x = element_text(margin = margin(t=10)),
    axis.title.y = element_text(margin = margin(r=10))
  )

ggsave(
  filename = here::here(release_dir, "hrs_age.png"),
  width = 14, height = 14, units = "cm"
)

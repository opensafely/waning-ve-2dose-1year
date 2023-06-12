library(tidyverse)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

tidy_dose3model <- read_csv(here::here("release20230605", "tidy_dose3model.csv")) 

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
      "Age*" = "age",
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
    filename = here::here("release20230605", glue::glue("{vars_group}.png")),
    width = 14, height = 20, units = "cm"
  )
  
  return(p)
  
}


create_plot("baseline")
create_plot("historyof")
create_plot("timeupdating")
  

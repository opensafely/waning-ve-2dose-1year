# setup ----
library(tidyverse)
library(glue)
library(survival)
library(splines)
library(broom)
library(broom.helpers)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# real list of model covariates
model_varlist <- readr::read_rds(
  here::here("analysis", "lib", "model_varlist.rds")
)

# create output directory
outdir <- here::here("output", "dose3", "model")
fs::dir_create(outdir)

# read data_timevarying
data_timevarying <- read_rds(here::here("output", "dose3", "data", "data_timevarying.rds")) 

# read baseline covariates
data_covs <- readr::read_rds(here::here("output", "data", "data_all.rds")) %>%
  right_join(data_timevarying %>% distinct(patient_id), by = "patient_id") %>%
  mutate(vax2_date = start_1_date - 14) %>%
  select(
    patient_id, 
    subgroup, jcvi_group, elig_date, vax2_date,
    vax2_brand = arm, 
    all_of(unname(unlist(model_varlist)))
  ) %>%
  mutate(across(vax2_brand, ~factor(.x)))

# print mean ages
data_covs %>%
  group_by(subgroup) %>%
  summarise(mean(age)) %>%
  print()

data_dose3model <- data_covs %>%
  left_join(data_timevarying, by = "patient_id") %>%
  # relevel imd so 3 is reference
  mutate(imd = fct_relevel(imd, levels(data_covs$imd)[3])) %>%
  mutate(across(subgroup, ~factor(.x, levels = subgroups))) %>%
  arrange(subgroup) %>%
  group_split(subgroup) %>%
  as.list()

data_dose3model %>%
  bind_rows() %>%
  mutate(across(subgroup, as.integer)) %>%
  group_by(subgroup, age) %>%
  count() %>%
  arrange(subgroup, age) %>%
  write_csv(file.path(outdir, "helperfile_agecounts.csv"))

if (FALSE) {
 # data summaries ----
for (s in seq_along(data_dose3model)) {
  
  cat(
    "\n--------------------------\n",
    "Subgroup: ",
    subgroups[s],
    "\n--------------------------\n"
  )
  
  cat("\nProportion of follow-up time spend in each status:\n")
  data_dose3model[[s]] %>%
    select(starts_with("status_")) %>%
    pivot_longer(cols = everything()) %>%
    group_by(name, value) %>%
    count() %>%
    ungroup() %>%
    pivot_wider(names_from = value, values_from = n, values_fill = 0) %>%
    mutate(prop1 = round(`1`/(`0`+`1`),3)) %>%
    print()
    
  cat("\nNumber of events in each status:\n")
  data_dose3model[[s]] %>%
    filter(ind_outcome==1) %>%
    select(starts_with("status_")) %>%
    pivot_longer(cols = everything()) %>%
    filter(value == 1) %>%
    group_by(name) %>% 
    count() %>% 
    print()
  
}

# define model formula ----
formula_dose3 <- formula(
  Surv(tstart, tstop, ind_outcome, type = "counting") ~ 
    # stratify by jcvi_group and elig_date, add vax2 date and age
    strata(jcvi_group, elig_date) + ns(vax2_date, 3) + poly(age, 2)
  ) %>%
  # add covariates
  update.formula(
    as.formula(
      str_c(
        ". ~ .",
        str_c(
          c(
            model_varlist$demographic[model_varlist$demographic != "age"],
            # - remove multimorb because including individual diagnoses
            # - remove pregnancy because only defined at baseline, 
            #   and difficult to accurately capture as time-varying variable,
            #   also only relevant for under 50s
            model_varlist$clinical[!(model_varlist$clinical %in% c("multimorb", "pregnancy"))],
            model_varlist$multimorb,
            "status_endoflife", "status_cancer", 
            "status_postest",
            "status_planned", "status_unplanned", "status_covidunplanned"
          ),
          collapse = "+"
          ),
        sep = "+"
        )
    )
  )

# print formula to check correct
cat("\n------\nModel formula:\n------\n")
print(formula_dose3)
cat("\nNote: `vax2_brand` added later for subgroups 1&2.\n")

# fit models separately within subgroups
glance <- list()
tidy <- list()
hrs_age <- list()
hrs_vax2_date <- list()
for (s in seq_along(data_dose3model)) {
  
  cat("\nSubgroup ", subgroups[s], ":\nstarted...\n")
  
  # add primary course brand to model if subgroups 1 or 2
  if (s %in% c(1,2)) {
    formula_subgroup <- formula_dose3 %>%
      update.formula(formula(. ~ . + vax2_brand))
  } else {
    formula_subgroup <- formula_dose3
  }
  
  # print warnings from coxph as they occur
  options(warn = 1)
  
  dose3model <- coxph(
    formula = formula_subgroup,
    data = data_dose3model[[s]],
    na.action = "na.fail"
  )
  
  write_rds(
    dose3model,
    file.path(outdir, glue("dose3model_{s}.rds")),
    compress = "gz"
  )
  
  # get predictions for vax2_date and age for plotting
  # code adapted from https://cran.r-project.org/web/packages/survival/vignettes/splines.pdf
  ptemp <- termplot(dose3model, se=TRUE, plot=FALSE, terms=1:2)
  
  # function for creating the data
  create_data <- function(var) {
    
    varterm <- as_tibble(ptemp[[var]])
    # get prediction for median age
    median_val <- median(data_dose3model[[s]][[var]])
    centre <- varterm %>% filter(x == median_val) %>% pull(y)
    data <- varterm %>%
      mutate(
        y_trans = y,
        conf.ll = y-(1.96*se),
        conf.ul = y+(1.96*se)
      ) %>%
      # rescale y relative to value at median_val
      mutate(across(c(y_trans, conf.ll, conf.ul), ~exp(.x - centre))) 
    
    attributes(data)$median_val <- median_val
    attributes(data)$var <- var
    
    return(data)
    
  }
  
  hrs_age[[s]] <- create_data("age")
  
  hrs_vax2_date[[s]] <- create_data("vax2_date")

  create_plot <- function(data) {
    data %>%
      ggplot(aes(x = x)) +
      geom_ribbon(aes(ymin = conf.ll, ymax = conf.ul), alpha = 0.2) +
      geom_line(aes(y = y_trans)) +
      geom_vline(xintercept = attributes(data)$median_val, linetype = "dashed") +
      scale_y_log10() +
      labs(x = attributes(data)$var, y = "HR") +
      theme_bw() 
  }
  
  hrs_age[[s]] %>% 
    create_plot() %>%
    ggsave(
      filename = file.path(outdir, glue("hrs_age_{s}.png"))
    )
  hrs_vax2_date[[s]] %>% 
    create_plot() %>%
    ggsave(
      filename = file.path(outdir, glue("hrs_vax2_date_{s}.png"))
    )
  
  # summary datasets
  glance[[s]] <-
    broom::glance(dose3model) %>%
    add_column(
      model = s,
      convergence = dose3model$info[["convergence"]],
      ram = format(object.size(dose3model), units="GB", standard="SI", digits=3L),
      .before = 1
    )
  
  tidy[[s]] <-
    broom.helpers::tidy_plus_plus(
      dose3model,
      exponentiate = TRUE,
    ) %>%
    transmute(
      model = s,
      variable, label, reference_row, 
      estimate, std.error, conf.low, conf.high, statistic, p.value
    )
  
  cat("...completed.\n")
  
}


# collapse outputs and save

glance %>%
  bind_rows() %>%
  write_csv(file.path(outdir, "glance_dose3model.csv"))

tidy %>%
  bind_rows() %>%
  write_csv(file.path(outdir, "tidy_dose3model.csv"))

hrs_age %>% 
  bind_rows(.id = "subgroup") %>%
  write_csv(
    file.path(outdir, "hrs_age.csv")
  )

hrs_vax2_date %>% 
  bind_rows(.id = "subgroup") %>%
  write_csv(
  file.path(outdir, "hrs_vax2_date.csv")
)
 
}
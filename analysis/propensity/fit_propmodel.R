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
outdir <- here::here("output", "propensity", "model")
fs::dir_create(outdir)

# read dataset
data_propmodel <- read_rds(here::here("output", "propensity", "data", "data_propmodel.rds")) 

# print mean ages
data_propmodel %>%
  group_by(subgroup) %>%
  summarise(mean(age)) %>%
  print()

data_propmodel <- data_propmodel %>%
  # age mean centred within subgroups
  group_by(subgroup) %>%
  mutate(age_smc = age - mean(age), .after = "age") %>%
  ungroup() %>%
  select(-age) %>%
  # relevel imd so 3 is reference
  mutate(imd = fct_relevel(imd, levels(data_propmodel$imd)[3])) %>%
  mutate(age_smc_squared = age_smc^2, .after = "age_smc") %>%
  group_split(subgroup) %>%
  as.list()

# define model formula ----

formula_prop <- formula(
  Surv(tstart, tstop, ind_outcome, type = "counting") ~ 
    # stratify by jcvi_group and elig_date
    strata(jcvi_group, elig_date) + ns(vax2_date, 3)
  ) %>%
  # add covariates
  update.formula(
    as.formula(
      str_c(
        ". ~ .",
        str_c(
          c(
            "age_smc", "age_smc_squared",
            model_varlist$demographic[model_varlist$demographic != "age"],
            # - remove multimorb because including individual diagnoses
            # - remove pregnancy because only defined at baseline, 
            #   and difficult to accurately capture as time-varying variable,
            #   also only relevant for under 50s
            model_varlist$clinical[!(model_varlist$clinical %in% c("multimorb", "pregnancy"))],
            model_varlist$multimorb,
            "endoflife", "postest",
            "inhosp_planned", "inhosp_unplanned", "inhosp_covidunplanned"
          ),
          collapse = "+"
          ),
        sep = "+"
        )
    )
  )

# print formula to check correct
print(formula_prop)

# fit models separately within subgroups
glance <- list()
tidy <- list()
for (s in subgroup_labels) {
  
  cat("Subgroup ", subgroups[s], ": started.\n")
  
  # add primary course brand to model if subgroups 1 or 2
  if (s %in% c(1,2)) {
    formula_subgroup <- formula_prop %>%
      update.formula(formula(. ~ . + vax2_brand))
  } else {
    formula_subgroup <- formula_prop
  }
  
  propmodel <- coxph(
    formula = formula_subgroup,
    data = data_propmodel[[s]],
    na.action = "na.fail"
  )
  
  write_rds(
    propmodel,
    file.path(outdir, glue("propmodel_{s}.rds")),
    compress = "gz"
  )
  
  glance[[s]] <-
    broom::glance(propmodel) %>%
    add_column(
      model = s,
      convergence = propmodel$info[["convergence"]],
      ram = format(object.size(propmodel), units="GB", standard="SI", digits=3L),
      .before = 1
    )
  
  tidy[[s]] <-
    broom.helpers::tidy_plus_plus(
      propmodel,
      exponentiate = TRUE
    ) %>%
    transmute(
      model = s,
      variable, label, reference_row, n_obs, n_event,
      estimate, std.error, conf.low, conf.high, statistic, p.value
    )
  
  cat("Subgroup ", subgroups[s], ": completed.\n")
  
}

glance %>%
  bind_rows() %>%
  write_csv(file.path(outdir, "glance_propmodel.csv"))

tidy %>%
  bind_rows() %>%
  write_csv(file.path(outdir, "tidy_propmodel.csv"))

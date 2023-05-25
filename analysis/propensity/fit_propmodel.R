# setup ----
library(tidyverse)
library(lubridate)
library(glue)
library(survival)
library(splines)


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

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# real list of model covariates
model_varlist <- readr::read_rds(
  here::here("analysis", "lib", "model_varlist.rds")
)

# output directory
outdir <- here::here("output", "propensity", "model")
fs::dir_create(outdir)

# define model formula ----

formula_prop <- formula(
  Surv(tstart, tstop, ind_outcome, type = "counting") ~ 
    # stratify by jcvi_group, elig_date and region
    strata(jcvi_group, elig_date, region) + ns(vax2_date, 3)
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
            model_varlist$clinical[model_varlist$clinical != "multimorb"],
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
print(formula_prop)
# Surv(tstart, tstop, ind_outcome, type = "counting") ~ strata(jcvi_group, 
#     elig_date, region) + ns(vax2_date, 3) + age_smc + age_smc_squared + 
#     sex + imd + ethnicity + bmi + learndis + sev_mental + flu_vaccine + 
#     test_hist_n + pregnancy + crd + chd + cld + ckd + cns + diabetes + 
#     immunosuppressed + endoflife + postest + inhosp_planned + 
#     inhosp_unplanned + inhosp_covidunplanned


# fit models separately within subgroups
for (s in subgroup_labels) {
  
  cat("Subgroup: ", subgroups[s], "\n")
  
  res_cox <- coxph(
    formula = formula_prop,
    data = data_propmodel[[s]],
    na.action = "na.fail"
  )
  
  print(summary(res_cox))
  
  write_rds(
    res_cox,
    file.path(outdir, glue("res_cox_{s}.rds")),
    compress = "gz"
  )
  
}

library(tidyverse)
library(yaml)
library(here)
library(glue)

# create folder for metadata
fs::dir_create(here::here("analysis", "lib"))

################################################################################
K <- 12L # the number of comparison periods

study_parameters <-
  list(
    seed = 123456L,
    n = 100000L, # number of individuals in dummy data
    K = K, 
    ref_age_1 = "2021-03-31", # reference date for calculating age for phase 1 groups
    ref_age_2 = "2021-07-01", # reference date for calculating age for phase 2 groups
    ref_cev = "2021-01-18", # reference date for calculating eligibility for phase 1 group 4 (CEV)
    ref_ar = "2021-02-15", # reference date for calculating eligibility for phase 1 group 5 (at-risk)
    pandemic_start = "2020-01-01", # rough start date for pandemic in UK
    start_date = "2020-12-08", # start of phase 1 vaccinations
    # start_date_pfizer = "2020-12-08",
    # start_date_az = "2021-01-04",
    # start_date_moderna = "2021-03-04",
    end_date = "2022-06-07", # day the study definition was run
    end_date_model = "2022-03-31" # based on availability of hospital data
  ) 

readr::write_rds(study_parameters, here::here("analysis", "lib", "study_parameters.rds"))
jsonlite::write_json(study_parameters, path = here::here("analysis", "lib", "study_parameters.json"), auto_unbox = TRUE, pretty=TRUE)

################################################################################
# jcvi_groups ----
jcvi_groups <- 
  tribble(
    ~group, ~definition,
    "01", "longres_group AND age_1 > 65",
    "02", "age_1 >=80",
    "03", "age_1 >=75",
    "04a", "age_1 >=70",
    "04b", "cev_group AND age_1 >=16",
    "05", "age_1 >=65",
    "06", "atrisk_group AND age_1 >=16",
    "07", "age_1 >=60",
    "08", "age_1 >=55",
    "09", "age_1 >=50",
    "10", "age_2 >=40",
    "11", "age_2 >=30",
    "12", "age_2 >=18",
    "99", "DEFAULT",
  )

readr::write_csv(jcvi_groups, here::here("analysis", "lib", "jcvi_groups.csv"))

################################################################################
# create elig_dates ----
# group elig_date if within 7 days of previous elig_date (within jcvi_group)
elig_dates <-
  tribble(
    ~date, ~description, ~jcvi_groups,
    "2020-12-08", "jcvi_group='01' OR jcvi_group='02'", "01, 02", 
    "2021-01-18", "jcvi_group='03' OR jcvi_group='04a' OR jcvi_group='04b'", "03, 04a, 04b",
    ###
    "2021-02-15", "jcvi_group='05' OR jcvi_group='06'", "05, 06",
    ###
    "2021-02-22", "age_1 >= 64 AND age_1 < 65", "07", 
    "2021-03-01", "age_1 >= 60 AND age_1 < 64", "07",
    ###
    # combine 2 rows as < 7 days between elig_dates
    "2021-03-08", "age_1 >= 55 AND age_1 < 60", "08",
    # "2021-03-08", "age_1 >= 56 AND age_1 < 60", "08",
    # "2021-03-09", "age_1 >= 55 AND age_1 < 56", "08",
    ###
    "2021-03-19", "age_1 >= 50 AND age_1 < 55", "09",
    ###
    "2021-04-13", "age_2 >= 45 AND age_1 < 50", "10",
    # combine 3 rows as < 7 days between elig_dates
    "2021-04-26", "age_2 >= 40 AND age_1 < 45", "10",
    # "2021-04-26", "age_2 >= 44 AND age_1 < 45", "10",
    # "2021-04-27", "age_2 >= 42 AND age_1 < 44", "10",
    # "2021-04-30", "age_2 >= 40 AND age_1 < 42", "10",
    ###
    # combine 2 rows as < 7 days between elig_dates
    "2021-05-13", "age_2 >= 36 AND age_2 < 40", "11",
    # "2021-05-13", "age_2 >= 38 AND age_2 < 40", "11",
    # "2021-05-19", "age_2 >= 36 AND age_2 < 38", "11",
    # combine 3 rows as < 7 days between elig_dates
    "2021-05-21", "age_2 >= 30 AND age_2 < 36", "11",
    # "2021-05-21", "age_2 >= 34 AND age_2 < 36", "11",
    # "2021-05-25", "age_2 >= 32 AND age_2 < 34", "11",
    # "2021-05-26", "age_2 >= 30 AND age_2 < 32", "11",
    ###
    "2021-06-08", "age_2 >= 25 AND age_2 < 30", "12",
    # combine 3 rows as < 7 days between elig_dates
    "2021-06-15", "age_2 >= 18 AND age_2 < 25", "12",
    # "2021-06-15", "age_2 >= 23 AND age_2 < 25", "12",
    # "2021-06-16", "age_2 >= 21 AND age_2 < 23", "12",
    # "2021-06-18", "age_2 >= 18 AND age_2 < 21", "12",
    "2100-12-31", "DEFAULT", "NA",
  ) 

readr::write_csv(elig_dates, here::here("analysis", "lib", "elig_dates.csv"))

################################################################################
# create regions ----
regions <- tribble(
  ~region, ~ratio,
  "North East", 0.1,
  "North West", 0.1,
  "Yorkshire and The Humber", 0.1,
  "East Midlands", 0.1,
  "West Midlands", 0.1,
  "East", 0.1,
  "London", 0.2,
  "South West", 0.1,
  "South East", 0.1
)

readr::write_csv(regions, here::here("analysis", "lib", "regions.csv"))

################################################################################
# varlists

# clinical variables for model
clinical <- c(
  "BMI" = "bmi",
  "Learning disability" = "learndis",
  "Serious mental illness" = "sev_mental",
  "Morbidity count" = "multimorb",
  "Flu vaccine in previous 5 years" = "flu_vaccine",
  "Number of SARS-CoV-2 tests between 2020-05-18 and min_elig_date" = "test_hist_n",
  "Pregnancy" = "pregnancy"
)

# extra clinical variables for summary table 
# (those used to define morbidity count)
multimorb <-c(
  "Chronic respiratory disease" = "crd", 
  "Chronic heart disease" = "chd", 
  "Chronic liver disease" = "cld", 
  "Chronic kidney disease" = "ckd", 
  "Chronic neurological disease" = "cns",
  "Diabetes" = "diabetes",
  "Immunosuppression" = "immunosuppressed"
)

demographic <- c(
  "Age" = "age",
  "Sex" = "sex",
  "IMD" = "imd",
  "Ethnicity" = "ethnicity"
)

readr::write_rds(
  list(demographic = demographic, clinical = clinical, multimorb = multimorb),
  here::here("analysis", "lib", "model_varlist.rds")
)

################################################################################
# strata vars for cox model ----
strata_vars <- c(
  "Region" = "region",
  "JCVI group" = "jcvi_group",
  "Date of eligibility for 1st dose" = "elig_date"
)

readr::write_rds(
  strata_vars,
  here::here("analysis", "lib", "strata_vars.rds")
)

################################################################################
# subgroups ----
subgroups <- c(
  "65+ years",
  "18-64 years and clinically vulnerable", 
  "40-64 years",
  "18-39 years"
  )

readr::write_rds(
  subgroups,
  here::here("analysis", "lib", "subgroups.rds")
)

################################################################################
# all outcomes:
outcomes <- c(
  "Any SARS-CoV-2 test" = "anytest", 
  "Positive SARS-CoV-2 test" = "postest", 
  "COVID-19 hospitalisation (APCS)" = "covidadmitted",
  "COVID-19 hospitalisation (ECDS)" = "covidemergency",
  "COVID-19 death" = "coviddeath", 
  "Non-COVID-19 death" = "noncoviddeath")

readr::write_rds(
  outcomes,
  here::here("analysis", "lib", "outcomes.rds")
)

outcomes_model <- outcomes
readr::write_rds(
  outcomes,
  here::here("analysis", "lib", "outcomes_model.rds")
)

outcomes <- outcomes[-4]
outcomes_model <- outcomes

################################################################################
comparisons <- c("BNT162b2", "ChAdOx1", "both")
readr::write_rds(
  comparisons,
  here::here("analysis", "lib", "comparisons.rds")
)

comparisons <- comparisons[-3]

# ################################################################################
# # create bash script for generating study definitions from template
# create_study_definitions <- 
#   str_c("# Run this script to create a study_definition for each k",
#         "",
#         str_c("for i in {1..",K,"}; do"),
#         "sed -e \"s;%placeholder_k%;$i;g\" ./analysis/study_definition_k.py > ./analysis/study_definition_$i.py;",
#         "done;", sep = "\n")
# 
# create_study_definitions %>%
#   writeLines(here::here("analysis/create_study_definitions.sh"))
# 
# # create study definitions from template_study_definition.py
# check_create <- try(processx::run(command="bash", args= "analysis/create_study_definitions.sh"))
# 
# if (inherits(check_create, "try-error")) stop("Study definitions not created.")


################################################################################
# create action functions ----

## generic action function ----
action <- function(
  name,
  run,
  dummy_data_file=NULL,
  arguments=NULL,
  needs=NULL,
  highly_sensitive=NULL,
  moderately_sensitive=NULL
){
  
  outputs <- list(
    highly_sensitive = highly_sensitive,
    moderately_sensitive = moderately_sensitive
  )
  outputs[sapply(outputs, is.null)] <- NULL
  
  action <- list(
    run = paste(c(run, arguments), collapse=" "),
    dummy_data_file = dummy_data_file,
    needs = needs,
    outputs = outputs
  )
  action[sapply(action, is.null)] <- NULL
  
  action_list <- list(name = action)
  names(action_list) <- name
  
  action_list
}


## create comment function ----
comment <- function(...){
  list_comments <- list(...)
  comments <- map(list_comments, ~paste0("## ", ., " ##"))
  comments
}


## create function to convert comment "actions" in a yaml string into proper comments
convert_comment_actions <-function(yaml.txt){
  yaml.txt %>%
    str_replace_all("\\\n(\\s*)\\'\\'\\:(\\s*)\\'", "\n\\1")  %>%
    #str_replace_all("\\\n(\\s*)\\'", "\n\\1") %>%
    str_replace_all("([^\\'])\\\n(\\s*)\\#\\#", "\\1\n\n\\2\\#\\#") %>%
    str_replace_all("\\#\\#\\'\\\n", "\n")
}

## actions that extract and process data ----
apply_model_fun <- function(
  subgroup_label,
  comparison,
  outcome
) {
  
  subgroup <- subgroups[subgroup_label]
  
  splice(
    comment(glue("{comparison}; {subgroup_label}; {outcome}")),
    comment("preflight checks"),
    action(
      name = glue("preflight_{comparison}_{subgroup_label}_{outcome}"),
      run = "r:latest analysis/comparisons/preflight.R",
      arguments = c(comparison, subgroup_label, outcome),
      needs = list(
        "data_covariates_process", 
        glue("data_tte_process_{comparison}")
        ),
      highly_sensitive = list(
        model_input = glue("output/preflight/data/model_input_{comparison}_{subgroup_label}_{outcome}*.rds")
      ),
      moderately_sensitive = list(
        # eventcheck_table = glue("output/preflight/tables/eventcheck_{comparison}_{subgroup_label}_{outcome}*.html"),
        preflight_report = glue("output/preflight/tables/preflight_report_{comparison}_{subgroup_label}_{outcome}*.txt")
      )
    ),
    comment("apply cox model"),
    action(
      name = glue("apply_model_cox_{comparison}_{subgroup_label}_{outcome}"),
      run = "r:latest analysis/comparisons/apply_model_cox_update.R",
      arguments = c(comparison, subgroup_label, outcome),
      needs = list(
        glue("preflight_{comparison}_{subgroup_label}_{outcome}")),
      highly_sensitive = list(
        modelnumber = glue("output/models_cox/data/model*_{comparison}_{subgroup_label}_{outcome}*.rds"),
        model_tidy_rds = glue("output/models_cox/data/modelcox_tidy_{comparison}_{subgroup_label}_{outcome}*.rds"),
        model_glance_rds = glue("output/models_cox/data/modelcox_glance_{comparison}_{subgroup_label}_{outcome}*.rds")
      ),
      moderately_sensitive = list(
        model_tidy_txt = glue("output/models_cox/temp/modelcox_tidy_{comparison}_{subgroup_label}_{outcome}.txt"),
        model_glance_txt = glue("output/models_cox/temp/modelcox_glance_{comparison}_{subgroup_label}_{outcome}.txt")
      )
    )
    
  )
}

# specify project ----

## defaults ----
defaults_list <- list(
  version = "3.0",
  expectations= list(population_size=100000L)
)

subgroup_labels <- seq_along(subgroups)

## actions ----
actions_list <- splice(
  
  comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          "DO NOT EDIT project.yaml",
          "This file is created by create_project.R",
          "Edit and run create_project.R to update the project.yaml",
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"
  ),
  
  comment("####################################", 
          "study definition",
          "####################################"),
  comment("generate dummy data for study_definition"),
  action(
    name = "dummy_data_vax",
    run = "r:latest analysis/dummy_data_vax.R",
    moderately_sensitive = list(
      dummy_data = "analysis/dummy_data_vax.feather"
    )
  ),
  
  comment("study definition"),
  action(
    name = "generate_study_population",
    run = "cohortextractor:latest generate_cohort --study-definition study_definition_vax --output-format feather",
    dummy_data_file = "analysis/dummy_data_vax.feather",
    needs = list("dummy_data_vax"),
    highly_sensitive = list(
      cohort = "output/input_vax.feather"
    )
  ),
  
  comment("####################################", 
          "preprocessing",
          "####################################"),
  
  comment("process data from study_definition"),
  action(
    name = "data_input_process",
    run = "r:latest analysis/preprocess/data_input_process.R",
    needs = list("dummy_data_vax", "generate_study_population"),
    highly_sensitive = list(
      data_wide_vax_dates = "output/data/data_wide_vax_dates.rds",
      data_processed = "output/data/data_processed.rds"
    ),
    moderately_sensitive = list(
      data_properties = "output/tables/data_*_tabulate.txt"
    )
  ),
  
  comment("apply eligiblity criteria from boxes a and b"),
  action(
    name = "data_eligible_ab",
    run = "r:latest analysis/preprocess/data_eligible_ab.R",
    needs = list("data_input_process"),
    highly_sensitive = list(
      data_eligible_a = "output/data/data_eligible_a.rds",
      data_eligible_b = "output/data/data_eligible_b.rds"
    ),
    moderately_sensitive = list(
      eligibility_count_ab = "output/tables/eligibility_count_ab.csv",
      group_age_ranges = "output/lib/group_age_ranges.csv"
    )
  ),
  
  comment("####################################", 
          "second_vax_period",
          "####################################"),
  comment("identify second vaccination time periods"),
  comment("create dataset for identifying second vaccination time periods"),
  action(
    name = "data_2nd_vax_dates",
    run = "r:latest analysis/second_vax_period/data_2nd_vax_dates.R",
    needs = list("data_input_process", "data_eligible_ab"),
    highly_sensitive = list(
      data_vax_plot = "output/second_vax_period/data/data_vax_plot.rds",
      second_vax_period_dates_rds = "output/second_vax_period/data/second_vax_period_dates.rds"
    ),
    moderately_sensitive = list(
      second_vax_period_dates_txt = "output/second_vax_period/tables/second_vax_period_dates.txt"
    )
  ),
  
  comment("plot second vaccination time periods"),
  action(
    name = "plot_2nd_vax_dates",
    run = "r:latest analysis/second_vax_period/plot_2nd_vax_dates.R",
    needs = list("data_eligible_ab", "data_2nd_vax_dates"),
    moderately_sensitive = list(
      plots_by_region = "output/second_vax_period/images/plot_by_region_*.png",
      plots_by_region_data = "output/second_vax_period/images/plot_by_region_*.txt"
    )
  ),
  
  comment("apply eligiblity criteria from boxes c, d and e"),
  action(
    name = "data_eligible_cde",
    run = "r:latest analysis/preprocess/data_eligible_cde.R",
    needs = list("data_input_process", "data_eligible_ab", "data_2nd_vax_dates"),
    highly_sensitive = list(
      data_eligible_e_vax = "output/data/data_eligible_e_vax.rds",
      data_eligible_e_unvax = "output/data/data_eligible_e_unvax.rds",
      data_eligible_e = "output/data/data_eligible_e.csv"
    ),
    moderately_sensitive = list(
      eligibility_count_cde = "output/tables/eligibility_count_cde.csv"
    )
  ),
  
  comment("####################################", 
          "study definition covs",
          "####################################"),
  
  action(
    name = "generate_covs_data",
    run = "cohortextractor:latest generate_cohort --study-definition study_definition_covs --output-format feather",
    needs = list("data_eligible_cde"),
    highly_sensitive = list(
      cohort = "output/input_covs.feather"
    )
  ),
  
  comment("####################################", 
          "process covariates data",
          "####################################"),
  
  comment("(includes anytest_date)"),
  action(
    name = "data_covariates_process",
    run = "r:latest analysis/preprocess/data_covariates_process.R",
    needs = splice(
      "data_input_process", 
      "data_eligible_cde", 
      "generate_covs_data"
      ),
    moderately_sensitive = list(
      data_min_max_fu_csv = "output/lib/data_min_max_fu.csv"
    ),
    highly_sensitive = list(
      data_covariates = "output/data/data_all.rds"
    )
  ),
  
  # comment("min and max follow-up dates for plots"),
  # action(
  #   name = "data_min_max_fu",
  #   run = "r:latest analysis/comparisons/data_min_max_fu.R",
  #   needs = list("data_covariates_process"),
  #   moderately_sensitive = list(
  #     data_min_max_fu_csv = "output/lib/data_min_max_fu.csv"
  #   )
  # ),
  
  comment("####################################",
          "subsequent vaccination", 
          "####################################"),
  
  comment("plot cumulative incidence of subsequent vaccination"),
  action(
    name = "plot_cumulative_incidence",
    run = "r:latest analysis/subsequent_vax/plot_cumulative_incidence.R",
    needs = list("data_covariates_process"),
    moderately_sensitive = list(
      ci_vax = "output/subsequent_vax/images/ci_vax.png",
      survtable = "output/subsequent_vax/tables/survtable_redacted.csv"
    )
  ),
  
  comment("####################################",
          "table 1 for report", 
          "####################################"),
  comment("create table 1 for all and for each subgroup"),
  action(
    name = "table1",
    run = "r:latest analysis/report/table1.R",
    needs = list(
      "data_eligible_ab",
      "data_eligible_cde",
      "data_covariates_process"
      ),
    moderately_sensitive = list(
      eligibility_count_p1 = "output/tables/eligibility_count_p1.csv",
      eligibility_count_all = "output/tables/eligibility_count_all.csv",
      table_csv = "output/report/tables/table1_*_REDACTED.csv",
      table_html = "output/report/tables/table1_*_REDACTED.html"
    )
  ),
  
  comment("####################################",
          "process time to event data", 
          "####################################"),
  
  comment(glue("process tte data")),
  splice(unlist(lapply(
    comparisons,
    function(x)
      action(
        name = glue("data_tte_process_{x}"),
        run = "r:latest analysis/comparisons/data_tte_process.R",
        arguments = x,
        needs = list("data_covariates_process"),
        highly_sensitive = list(
          data_tte_brand_outcome = glue("output/tte/data/data_tte_{x}*.rds")
        ),
        moderately_sensitive = list(
          event_counts_csv = glue("output/tte/data/event_counts_{x}.csv"),
          event_counts_txt = glue("output/tte/tables/event_counts_{x}.txt")
        )
      )
  ), recursive = FALSE)),
  

  comment("check distribution of follow-up time in relation to variant dates"),
  splice(
    unlist(
      lapply(
        seq_along(subgroups),
        function(x)
          action(
            name = glue("check_fu_{x}"),
            run = "r:latest analysis/comparisons/check_fu.R",
            arguments = x,
            needs = list(
              "data_covariates_process", 
              "data_tte_process_BNT162b2", 
              "data_tte_process_ChAdOx1"
              ),
            moderately_sensitive = list(
              check_fu_plot = glue("output/tte/images/check_fu_{x}.png"),
              check_fu_plot_data = glue("output/tte/images/check_fu_{x}.csv")
            )
          )
      ),
      recursive = FALSE
    )
  ),
  
  comment("####################################",
          "apply models", 
          "####################################"),
  splice(
    # over subgroups
    unlist(lapply(
      comparisons,

      function(x) {
        if (!(x %in% "BNT162b2")) {
          ys <- subgroup_labels[subgroups != "18-39 years"]
        } else {
          ys <- subgroup_labels
        }
        unlist(lapply(
          ys,
          function(y)
            splice(
            unlist(lapply(
              unname(outcomes_model),

              function(z)
              apply_model_fun(
                comparison = x,
                subgroup_label = y,
                outcome = z)
            ),
            recursive = FALSE)
            )
          
        ),
        recursive = FALSE)
      }
    ), recursive = FALSE)
  ),
  
  comment("combine all estimates for release"),
  action(
    name = glue("combine_estimates"),
    run = "r:latest analysis/comparisons/combine_estimates.R",
    needs = splice(
      lapply(
        comparisons[comparisons != "both"], 
        function(x) glue("data_tte_process_{x}")
        ),
      as.list(unlist(lapply(
        comparisons,
        function(x)
        {
          if (x %in% c("ChAdOx1", "both")) {
            ys <- subgroup_labels[subgroups != "18-39 years"]
          } else {
            ys <- subgroup_labels
          }
          unlist(lapply(
            ys,
            function(y)
              unlist(lapply(
                unname(outcomes_model),
                function(z)
                  glue("apply_model_cox_{x}_{y}_{z}")
              ), recursive = FALSE)
          ), recursive = FALSE)
        }
      ), recursive = FALSE))),
    moderately_sensitive = list(
      event_counts_all = glue("output/release_objects/event_counts*.csv"),
      estimates_all = glue("output/release_objects/estimates*.csv")
    )
  ),
  
  comment("####################################",
          "tables for appendix", 
          "####################################"),
  action(
    name = "appendix_table",
    run = "r:latest analysis/post_release/appendix_table.R",
    needs = list(
      "combine_estimates"
      ),
    moderately_sensitive = list(
      plot_check = "output/release_objects/appendix_table*.csv"
    )
  ),
  
  comment("####################################",
          "plot for manuscript", 
          "####################################"),
  action(
    name = "plot_cox_all",
    run = "r:latest analysis/post_release/plot_cox_all.R",
    needs = splice(
      "combine_estimates",
      "plot_cumulative_incidence",
      lapply(1:4, function(x) glue("check_fu_{x}"))
    ),
    moderately_sensitive = list(
      hr_vax_ci = "output/release_objects/images/hr_vax_ci.png"
    )
  )
  
)

## combine everything ----
project_list <- splice(
  defaults_list,
  list(actions = actions_list)
)

## convert list to yaml, reformat comments and whitespace,and output ----
as.yaml(project_list, indent=2) %>%
  # convert comment actions to comments
  convert_comment_actions() %>%
  # add one blank line before level 1 and level 2 keys
  str_replace_all("\\\n(\\w)", "\n\n\\1") %>%
  str_replace_all("\\\n\\s\\s(\\w)", "\n\n  \\1") %>%
  writeLines(here("project.yaml"))


## grab all action names and send to a txt file

names(actions_list) %>% tibble(action=.) %>%
  mutate(
    model = str_detect(action, "model"),
    model_number = cumsum(model)
  ) %>%
  group_by(model_number) %>%
  summarise(
    sets = paste(action, collapse=" ")
  ) %>% pull(sets) %>%
  paste(collapse="\n") %>%
  writeLines(here("actions.txt"))

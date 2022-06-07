library(tidyverse)
library(glue)
library(lubridate)
library(gt)

################################################################################
fs::dir_create(here::here("output", "tables"))
fs::dir_create(here::here("output", "report", "tables"))

################################################################################
# read list of covariates for model
model_varlist <- readr::read_rds(
  here::here("analysis", "lib", "model_varlist.rds"))

# read strata_vars
strata_vars <- readr::read_rds(
  here::here("analysis", "lib", "strata_vars.rds"))
strata_vars <- strata_vars[strata_vars!="elig_date"]

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

################################################################################
# processed data
data_all <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) 

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

# function to be applied in dplyr::filter
no_evidence_of <- function(cov_date, index_date) {
  is.na(cov_date) | index_date < cov_date
}

################################################################################
# prepare data
data_tables_0 <- data_all %>%
  mutate(group = if_else(arm == "unvax", "unvax", "vax"))

# eligible for comparison period 1
eligibility_count <- data_tables_0 %>% 
  group_by(group) %>%
  count() %>%
  ungroup() %>%
  transmute(
    description = glue("{group}: satisfying eligibility criteria up to and including box E."),
    n
  )

# remove if death before start of comparison 1 (includes day 0)
data_tables <- data_tables_0 %>%
  filter_at(
    vars("death_date"),
    all_vars(no_evidence_of(., start_1_date))) 

eligibility_count <- eligibility_count %>%
  bind_rows(
    data_tables %>% 
      group_by(group) %>%
      count() %>%
      ungroup() %>%
      transmute(
        description = glue("{group}: after removing those who died before start of period 1."),
        n
      ))

# remove if dereg before start of comparison 1
data_tables <- data_tables %>%
  filter_at(
    vars("dereg_date"),
    all_vars(no_evidence_of(., start_1_date))) 

eligibility_count <- eligibility_count %>%
  bind_rows(
    data_tables %>% 
      group_by(group) %>%
      count() %>%
      ungroup() %>%
      transmute(
        description = glue("{group}: after removing those who deregistered before start of period 1."),
        n
      ))

# remove if subsequent_vax before start of comparison 1
data_tables <- data_tables %>%
  filter_at(
    vars("subsequent_vax_date"),
    all_vars(no_evidence_of(., start_1_date))) 

eligibility_count <- eligibility_count %>%
  bind_rows(
    data_tables %>% 
      group_by(group) %>%
      count() %>%
      ungroup() %>%
      transmute(
        description = glue("{group}: after removing those who received a subsequent dose before start of period 1."),
        n
      ))

eligibility_count_p1 <- eligibility_count %>%
  mutate(group = str_extract(description, "\\w+:")) %>%
  arrange(group) %>%
  # round up to nearest 7
  mutate(across(n, ~ceiling_any(.x, to=7))) %>%
  group_by(group) %>%
  mutate(n_removed = lag(n) - n) %>%
  ungroup() %>%
  select(-group)

readr::write_csv(
  eligibility_count_p1,
  here::here("output", "tables", "eligibility_count_p1.csv"))

################################################################################
# combine eligibility count tables
eligibility_count_all <- bind_rows(
  readr::read_csv(here::here("output", "tables", "eligibility_count_ab.csv")),
  readr::read_csv(here::here("output", "tables", "eligibility_count_cde.csv")),
  eligibility_count_p1 %>% mutate(stage="p1")
)

readr::write_csv(
  eligibility_count_all,
  here::here("output", "tables", "eligibility_count_all.csv"))

################################################################################
# split data in subgroups
data_tables <- data_tables %>%
  select(patient_id, arm, region, jcvi_group, subgroup,
         all_of(unname(unlist(model_varlist)))) %>% 
  group_split(subgroup)

################################################################################
# function to summarise each variable and (within variables) 
# round all frequencies up to nearest 7
summary_var <- function(.data, var) {
  out <- .data %>%
    group_by(arm, !! sym(var)) %>%
    count() %>%
    ungroup(!! sym(var)) %>%
    mutate(arm_total = sum(n)) %>%
    ungroup() %>%
    mutate(across(n, ~ceiling_any(.x, to=7))) %>%
    mutate(percent = round(100*n/arm_total,0)) %>%
    mutate(across(n, ~scales::comma(.x, accuracy = 1))) %>%
    mutate(value = as.character(glue("{n} ({percent}%)"))) %>%
    select(arm, !! sym(var), value) %>%
    pivot_wider(
      names_from = arm, 
      values_from = value
    ) %>%
    mutate(variable = var) %>%
    rename("category" = var)
  
  if (is.logical(out$category)) {
    out <- out %>%
      filter(category) %>%
      mutate(across(category, ~ "yes"))
  }
  
  out %>% mutate(across(category, as.character))
  
}

################################################################################
# make table1 for all and each subgroup
for (i in c(seq_along(data_tables))) {
  cat(glue("---- loop {i} ----"), "\n")
  cat("---- define obejcts ----\n")
  variables <- c(unname(strata_vars), unname(unlist(model_varlist)))
  variables <- variables[variables != "age"]
  vars_ordered_levs <- c("region", "jcvi_group", "sex", "imd", "ethnicity", "bmi", "multimorb", "test_hist_n")
  
  # tibble for assigning tidy variable names
  var_labels <- tibble(
    variable = c(
      strata_vars,
      model_varlist$clinical, 
      model_varlist$multimorb, 
      model_varlist$demographic),
    variable_label = names(c(
      strata_vars, 
      model_varlist$clinical, 
      model_varlist$multimorb,
      model_varlist$demographic))
  )
  
  if (i == 0) {
    
    data <- bind_rows(data_tables) %>%
      droplevels()
    
    subgroup <- "All subgroups"
    subgroup_label <- 5
    
    variables <- c("subgroup", variables)
    vars_ordered_levs <- c("subgroup", vars_ordered_levs)
    
    var_labels <- var_labels %>%
      add_row(variable = "subgroup", variable_label = "Subgroup",
              .before=TRUE)
    
    min_elig_date <- "2020-12-08"
    
  } else {
    
    data <- data_tables[[i]] %>%
      droplevels()
    
    subgroup <- unique(data$subgroup)
    subgroup_label <- which(subgroups == subgroup)
    
    min_elig_date <- data_all %>%
      filter(subgroup %in% subgroup) %>%
      summarise(min_elig_date = min(elig_date))
    min_elig_date <- min_elig_date$min_elig_date
    
  }
  
  # function for creating tibble of categories for each variable
  var_tibble <- function(var) {
    var_region <- tibble(
      variable = var,
      category = levels(data[[var]])
    )
  }
  # tibble for specifying order of variables and categories
  var_order <- tibble(
    variable = variables
  ) %>%
    left_join(
      bind_rows(
        lapply(vars_ordered_levs,
               var_tibble)),
      by = "variable"
    ) %>%
    mutate(across(category, ~ if_else(is.na(.x), "yes", .x)))
  
  cat("---- make table 1 ----\n")
  # make table1
  table1 <- bind_rows(lapply(
    variables,
    function(x)
      data %>% summary_var(var = x)
  ))
  
  cat("---- tidy table 1 ----\n")
  # variables under "Evidence of" heading
  history_of_vars <- c(
    "crd",
    "chd", 
    "cld", 
    "ckd", 
    "cns",
    "diabetes",
    "immunosuppressed",
    # "asplenia",
    "learndis", 
    "sev_mental")
  # tidy table1
  table1_tidy <- var_order %>% 
    left_join(var_labels, by = "variable") %>%
    left_join(table1, by = c("category", "variable")) %>%
    mutate(across(category,
                  ~ if_else(variable %in% history_of_vars, variable_label, .x))) %>%
    mutate(across(variable_label,
                  ~ if_else(variable %in% history_of_vars, "Evidence of", .x))) %>%
    mutate(across(variable_label, ~ str_replace(.x, "min_elig_date", as.character(min_elig_date)))) %>%
    rename(Variable = variable_label, Characteristic = category, Unvaccinated = unvax) %>%
    select(-variable) %>%
    select(Variable, Characteristic, everything()) %>%
    mutate(across(c(BNT162b2, ChAdOx1, Unvaccinated), 
                  ~ if_else(is.na(.x), "- (-%)", .x))) 
  
  # age summary
  age_summary <- data %>%
    group_by(arm) %>%
    summarise(
      median = median(age, na.rm=TRUE),
      iqr = IQR(age, na.rm = TRUE),
      .groups = "keep") %>% 
    ungroup() %>%
    transmute(arm, value = as.character(glue("{median} ({iqr})"))) %>%
    pivot_wider(names_from = "arm", values_from = "value") %>%
    mutate(Variable = "Age", Characteristic = "Median (IQR)") 
  
  # age_missing <- data %>% 
  #   group_by(arm) %>%
  #   summarise(
  #     missing = sum(is.na(age)), 
  #     .groups = "keep") %>%
  #   ungroup() %>%
  #   transmute(arm, value = scales::comma(missing, accuracy = 1)) %>%
  #   pivot_wider(names_from = "arm", values_from = "value") %>%
  #   mutate(Variable = "Age", Characteristic = "Missing") 
  
  table1_tidy_n <- data %>% 
    group_by(arm) %>% 
    count() %>% 
    ungroup() %>% 
    mutate(across(n, ~ceiling_any(.x, to=7))) %>%
    pivot_wider(names_from = arm, values_from = n) %>% 
    mutate(Variable = "", Characteristic = "N") %>%
    mutate(across(c(BNT162b2, ChAdOx1, unvax), 
                  ~ scales::comma(.x, accuracy = 1))) %>%
    bind_rows(
      age_summary
      ) %>%
    rename(Unvaccinated = unvax) %>%
    bind_rows(
      table1_tidy
    )
  
  cat("---- save table1.csv ----\n")
  # save table1_tidy
  readr::write_csv(table1_tidy_n,
                   here::here("output", "report", "tables", glue("table1_{subgroup_label}_REDACTED.csv")))
  
  cat("---- save table1.html ----\n")
  table1_tidy_n %>%
    gt(
      groupname_col="Variable",
      rowname_col = "Characteristic"
    ) %>%
    tab_header(
      title = glue("Subgroup: {subgroup}"),
      subtitle = "Patient characteristics as of second vaccination period + 2 weeks") %>%
    tab_style(
      style = cell_text(weight="bold"),
      locations = list(
        cells_column_labels(
          columns = everything()
        ),
        cells_row_groups(
          groups = everything()
        ))
    ) %>%
    gtsave(
      filename = glue("table1_{subgroup_label}_REDACTED.html"),
      path = here::here("output", "report", "tables")
    )
  
}







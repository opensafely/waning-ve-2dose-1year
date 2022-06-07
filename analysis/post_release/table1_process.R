library(tidyverse)
library(glue)
library(flextable)
library(officer)
library(magrittr)
library(kableExtra)

################################################################################
if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")

################################################################################
# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)


################################################################################
# prepare the table 1 data
table_out0 <- bind_rows(
  lapply(
    1:4,
    function(x) 
      readr::read_csv(here::here(release_folder, "table1", glue("table1_{x}_REDACTED.csv")),
                      show_col_types = FALSE) %>%
      mutate(subgroup = x)
  ) 
) %>%
  pivot_wider(
    names_from = subgroup,
    values_from = c(BNT162b2, ChAdOx1, Unvaccinated),
    names_glue = "{subgroup}_{.value}"
  ) %>%
  select(-`4_ChAdOx1`) %>%
  select(Variable, Characteristic, starts_with(as.character(subgroup_labels))) %>%
  mutate(across(starts_with(as.character(subgroup_labels)),
                ~ if_else(
                  is.na(.x),
                  "-",
                  .x))) %>%
  mutate(across(Variable, ~if_else(Characteristic == "N", "N", .x))) %>%
  mutate(across(Characteristic, ~if_else(.x == "N", "", .x))) 

################################################################################
# prepare the column names
# column names for table
table_names <- names(table_out0)
# remove subgroup label
col_names <- str_remove(table_names, "\\d_")
names(col_names) <- table_names
# number of columns for each subgroup
cols_subtype <- sapply(
  as.character(subgroup_labels), 
  function(x) sum(str_detect(table_names, glue("{x}_")))
)
names(cols_subtype) <- subgroups
# reorder
cols_subtype <- cols_subtype[subgroup_labels]

################################################################################
# create and save the version for the manuscript
variable_order_manuscript <- c(
  "N" = "N", 
  "Age" = "Age", 
  "Sex" = "Sex", 
  "IMD" = "IMD (1 is most deprived)", 
  "Ethnicity" = "Ethnicity",
  "BMI" = "BMI", 
  "Morbidity count" = "Morbidity count",
  "Number of SARS-CoV-2 tests between 2020-05-18 and 2020-12-08" = "Number of SARS-CoV-2 tests", 
  "Flu vaccine in previous 5 years" = "Flu vaccine")

col_names[1] <- "Characteristic"

table1_data_manuscript <- tibble(
  variable_long = names(variable_order_manuscript),
  variable_short = unname(variable_order_manuscript)
  ) %>%
  left_join(table_out0, by = c("variable_long" = "Variable")) %>%
  select(-variable_long) %>%
  rename(Variable = variable_short) %>%
  mutate(across(Characteristic, ~str_remove(.x, "\\s\\w+\\sdeprived"))) %>%
  mutate(across(Characteristic, 
                ~if_else(
                  Variable == "BMI",
                  str_remove_all(str_extract(.x, "\\(\\d.+\\)"), "\\(|\\)"),
                  .x
                  ))) %>%
  mutate(across(Characteristic, 
                ~if_else(
                  Variable == "BMI" & is.na(.x),
                  "<30",
                  .x
                ))) 

# create table1_manuscript.docx
page_width_docx <- 26 #cm
cell_padding <- 0 # this is a guess, not sure what default is
col1_with <- 1.75
col2_width <- 1.5
coli_width <- (page_width_docx - cell_padding*ncol(table1_data_manuscript) - col1_with - col2_width)/(ncol(table1_data_manuscript) - 2)

flextable1 <- table1_data_manuscript %>%
  flextable() %>%
  set_header_labels(
    values = as.list(col_names)
  ) %>%
  merge_v(j=1, part = "body") %>%
  merge_at(i=1,j=1:2, part="header") %>%
  add_header_row(
    values = c("", names(cols_subtype)),
    colwidths = c(2, unname(cols_subtype))
  ) %>%
  width(j=1, width = col1_with, unit = "cm") %>%
  width(j=2, width = col1_with, unit = "cm") %>%
  width(j=3:ncol(table1_data_manuscript), width = coli_width, unit = "cm") %>%
  #TODO
  # footnote() 
  fontsize(size = 8, part = "all") %>%
  theme_booktabs() %>%
  padding(
    padding.top = 1,
    padding.bottom = 1,
    part = "all"
  )

doc <- read_docx() %>%
  body_add_flextable(value = flextable1, split = FALSE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print(target = here::here(release_folder, "my_table.docx"))

################################################################################
# create and save the version for the supplement

# add footnote marker to not clinically vulnerable subgroups for supplement
names(cols_subtype)[3:4] <- str_c(
  names(cols_subtype)[3:4],
  footnote_marker_alphabet(1, format = "latex", double_escape = TRUE)
)

# create the version for the appendix
variable_order_supplement <- c(
  "N",
  "Region", 
  "JCVI group", 
  "Evidence of", 
  "Pregnancy")

table1_data_supplement <- tibble(Variable = variable_order_supplement) %>%
  left_join(table_out0, by = "Variable") %>%
  mutate(tmp = row_number()) %>%
  group_by(Variable) %>%
  mutate(tmp = mean(tmp)) %>%
  arrange(Characteristic,.by_group = TRUE) %>%
  ungroup() %>%
  arrange(tmp) %>%
  select(-tmp) %>%
  mutate(across(Characteristic, 
                ~if_else(.x == "Immunosuppression",
                         "Immunosupp- ression",
                         .x))) %>%
  mutate(across(Variable, ~if_else(is.na(.x), " ", .x))) %>%
  mutate(across(starts_with(as.character(1:4)), ~str_replace(.x, "- \\(-\\%\\)", "-")))


table1_supplement <- list(
  col_names = col_names,
  cols_subtype = cols_subtype,
  data = table1_data_supplement
)

readr::write_rds(
  table1_supplement,
  here::here(release_folder, "table1_supplement.rds"))

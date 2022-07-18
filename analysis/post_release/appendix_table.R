library(tidyverse)
library(glue)

cat("create and define release folder")
fs::dir_create(here::here("output", "release_objects"))
if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))
outcomes <- outcomes[!(outcomes %in% c("covidemergency", "anytest"))]
old_names <- names(outcomes)
new_names <- str_replace(old_names, "Positive", "positive")
new_names <- str_replace(new_names, "Non", "non")
new_names <- str_remove(new_names, " \\(APCS\\)")
names(outcomes) <- new_names

outcomes_order <- c(2,3,1,4)

cat("read event counts data")
event_counts <- readr::read_csv(file.path(release_folder, "event_counts_all.csv")) %>%
  mutate(across(subgroup, as.integer)) %>%
  filter(!is.na(subgroup)) %>%
  select(-person_years) %>%
  pivot_wider(
    names_from = arm,
    values_from = c(n, events)
  )

cat("read estimates data")
estimates_k <- readr::read_csv(file.path(release_folder, "estimates_all.csv")) %>%
  filter(variable == "k", !reference_row) %>%
  mutate(across(c("estimate", "conf.low", "conf.high"),
                ~format(round(exp(.x), 2), nsmall=2))) %>%
  filter(comparison != "both") %>%
  transmute(
    subgroup, 
    arm = comparison,
    outcome, 
    k = period,
    model,
    value = glue("{estimate} ({conf.low},{conf.high})")
    ) %>%
  pivot_wider(
    names_from = model, 
    values_from = value
  ) %>%
  pivot_wider(
    names_from = arm,
    values_from = c(unadjusted, adjusted)
  )


data <- event_counts %>%
  left_join(
    estimates_k,
    by = c("subgroup", "outcome", "k")
  )

appendix_table_docx <- function(s,o) {
  
  subgroup_long <- subgroups[s]
  outcome_long <- names(outcomes)[outcomes==o]
  
  data_select <- data %>%
    filter(subgroup==s, outcome == o) 
  
  if (s %in% 1:2) {
    data_table <- data_select %>%
      select(
        k, n_unvax, events_unvax,
        n_BNT162b2, events_BNT162b2, unadjusted_BNT162b2, adjusted_BNT162b2,
        n_ChAdOx1, events_ChAdOx1, unadjusted_ChAdOx1, adjusted_ChAdOx1
      )
  } else if (s %in% 3) {
    data_table <- data_select %>%
      select(
        k, n_unvax, events_unvax,
        n_ChAdOx1, events_ChAdOx1, unadjusted_ChAdOx1, adjusted_ChAdOx1
      )
  } else if (s %in% 4) {
    data_table <- data_select %>%
      select(
        k, n_unvax, events_unvax,
        n_BNT162b2, events_BNT162b2, unadjusted_BNT162b2, adjusted_BNT162b2
      ) 
  }
  
  if(Sys.getenv("OPENSAFELY_BACKEND") %in% "") {
    
    names_data_table <- names(data_table)
    # first header row
    table_colnames_1 <- str_remove(names_data_table, "n_|events_|unadjusted_|adjusted_")
    table_colnames_1 <- str_replace(table_colnames_1, "unvax", "Unvaccinated")
    table_colnames_1[1] <- ""
    table_colnames_1_unique <- unique(table_colnames_1)
    table_colnames_1_length <- sapply(
      lapply(table_colnames_1_unique, function(x) table_colnames_1 == x),
      sum
    )
    
    names(table_colnames_1) <- names_data_table
    # second header row
    table_colnames_2 <- str_remove(names_data_table, "_unvax|_BNT162b2|_ChAdOx1")
    table_colnames_2 <- str_replace(table_colnames_2, "adjusted", "adjusted HR")
    names(table_colnames_2) <- names_data_table
    
    cell_width <- vector(mode = "integer", length = length(table_colnames_2))
    cell_width[table_colnames_2 %in% "k"] <- 2L
    cell_width[table_colnames_2 %in% c("n", "events")] <- 2L
    cell_width[table_colnames_2 %in% c("unadjusted HR", "adjusted HR")] <- 3L
    
    border <- officer::fp_border()
    border_j <- c("k", "events_unvax") 
    if (s %in% 1:2) border_j[3] <- "adjusted_BNT162b2"
    
    flextable_out <- data_table %>%
      flextable::flextable() %>%
      flextable::set_header_labels(
        values = as.list(table_colnames_2)
      ) %>%
      flextable::add_header_row(
        values = table_colnames_1_unique,
        colwidths = table_colnames_1_length
      ) %>%
      flextable::align(i=1, align="left",part="header") %>%
      flextable::width(j=1:ncol(data_table), width = cell_width, unit = "cm") %>%
      flextable::border_remove() %>%
      flextable::hline(i = 1:2, border = border, part = "header") %>%
      flextable::vline(j = border_j, border = border, part = "all") %>%
      flextable::fontsize(size = 8, part = "all") %>%
      # theme_booktabs() %>%
      flextable::padding(
        padding.top = 1,
        padding.bottom = 1,
        part = "all"
      ) %>%
      flextable::set_caption(
        caption = glue("Counts and hazard ratios (HR) for {outcome_long} the {subgroup_long} subgroup"),
        autonum = autonum
      )
    
    doc <- flextable::body_add_flextable(doc, value = flextable_out, split = FALSE)  %>%
      officer::body_add_par(value = " ", style = NULL, pos = "after") # add blank line
    
  } else {
    
    readr::write_csv(
      data_table,
      file.path(release_folder, glue("appendix_table_{s}_{o}.csv"))
    )
    
  }
  
}

if(Sys.getenv("OPENSAFELY_BACKEND") %in% "") {
  
  autonum <- officer::run_autonum(
    seq_id = "table",
    pre_label = "Supplementary Table ",
    post_label = ": ",
    bkm = NULL,
    bkm_all = FALSE,
    prop = NULL,
    start_at = NULL,
    tnd = 0,
    tns = "-"
  )
  
  doc <- officer::read_docx() 
}


cat("apply tables function")
for (s in 1:4) {
  for (o in outcomes[outcomes_order]) {
    
    appendix_table_docx(s,o)
    
  }
}

if(Sys.getenv("OPENSAFELY_BACKEND") %in% "") {
  doc <- officer::body_end_section_landscape(doc)
  doc <- print(doc, target = file.path(release_folder, "appendix_table.docx"))
}




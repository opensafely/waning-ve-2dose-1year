#################

# This function takes a dataset, summarises the variables using:
#   * skimr::skim(),
#   * class(),
#   * and ,
# and saves the output to a .txt file

#################

data_properties <- function(
  data, # input data
  path, # path for output file
  coltypes = FALSE, # save file that summarises coltypes?
  skim = FALSE, # save a file with summary stats from skim?
  tabulate = TRUE # save a file of redacted summary stats?
  ) {
  
  # preliminaries
  library('tidyverse')
  source(here::here("analysis", "functions", "redaction_functions.R"))
  
  # reset options when exit the function
  op <- options()
  on.exit(options(op))
  
  # Output summary .txt
  options(width=200) # set output width for capture.output
  
  # prefix for output files
  filenamebase <- deparse(substitute(data))
  
  if (skim) {
    ## High-level variable overview ----
    capture.output(
      skimr::skim_without_charts(data),
      file = file.path(path, paste0(filenamebase, "_skim", ".txt")),
      split=FALSE
    )
  }
  
  if (coltypes) {
    ## list of column types ----
    capture.output(
      lapply(data, class),
      file = file.path(path, paste0(filenamebase, "_coltypes", ".txt"))
    )
  }
  
  ## tabulated data ----
  if (tabulate) {
    # delete file if it exists
    if(file.exists(file.path(path, paste0(filenamebase, "_tabulate", ".txt")))){
      file.remove(file.path(path, paste0(filenamebase, "_tabulate", ".txt")))
    }
    
    ### categorical and logical ----
    sumtabs_cat <-
      data %>%
      select(-ends_with("_id")) %>%
      select(where(is.character), where(is.logical), where(is.factor)) %>%
      map(redacted_summary_cat) %>%
      enframe()
    
    capture.output(
      walk2(sumtabs_cat$value, sumtabs_cat$name, print_cat),
      file = file.path(path, paste0(filenamebase, "_tabulate", ".txt")),
      append=FALSE
    )
    
    ### numeric ----
    sumtabs_num <-
      data %>%
      select(-ends_with("_id")) %>%
      select(where(~ {!is.logical(.x) & is.numeric(.x) & !is.Date(.x)})) %>%
      map(redacted_summary_num) %>%
      enframe()
    
    capture.output(
      walk2(sumtabs_num$value, sumtabs_num$name, print_num),
      file = file.path(path, paste0(filenamebase, "_tabulate", ".txt")),
      append=TRUE
    )
    
    ### dates ----
    
    sumtabs_date <-
      data %>%
      select(-ends_with("_id")) %>%
      select(where(is.Date)) %>%
      map(redacted_summary_date) %>%
      enframe()
    
    capture.output(
      walk2(sumtabs_date$value, sumtabs_date$name, print_num),
      file = file.path(path, paste0(filenamebase, "_tabulate", ".txt")),
      append=TRUE
    )
  }
}


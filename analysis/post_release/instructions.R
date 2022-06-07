# define release folder
release_folder <- here::here("release20220505")

### Flow chart
# process data for flow chart
source(here::here("analysis", "post_release", "flow.R")) 
rm(list = ls()[ls() != "release_folder"])

### Table 1
# process data for summary table in manuscript and supplementary material
source(here::here("analysis", "post_release", "table1_process"))
rm(list = ls()[ls() != "release_folder"])
# render table for manuscript
rmarkdown::render(
  here::here("analysis","post_release", "table1_process.Rmd"),
  knit_root_dir = release_folder,
    output_file = here::here(release_folder, "table1_process.docx"))
rm(list = ls()[ls() != "release_folder"])

### Subsequent vaccination
# plot cumulative incidence of subsequent vaccination
source(here::here("analysis", "subsequent_vax", "plot_cumulative_incidence.R"))
rm(list = ls()[ls() != "release_folder"])

### Metaregression
# preprocess data for metaregression
source(here::here("analysis", "post_release", "data_metareg_process.R"))
rm(list = ls()[ls() != "release_folder"])
# Run waning_metareg.do in Stata
# postprocess metaregression results
source(here::here("analysis", "post_release", "data_metareg_k.R"))
rm(list = ls()[ls() != "release_folder"])

### Results plots
# all hazard ratio plots
source(here::here("analysis", "post_release", "plot_cox_all.R"))
rm(list = ls()[ls() != "release_folder"])

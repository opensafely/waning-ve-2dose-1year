# render tables

rmarkdown::render(
  here::here("analysis","post_release", "table1_process.Rmd"),
  output_file = here::here("release20220221", "table1_process.docx"))

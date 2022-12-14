library(tidyverse)
library(glue)

# 
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))


# define variant dates
delta_start <- as.Date("2021-06-01")
omicron_start <- as.Date("2021-12-01")
delta_end <- as.Date("2021-12-15")

# labels for k
K <- study_parameters$K
ends <- seq(2, (K+1)*4, 4)
starts <- ends + 1
weeks_since_2nd_vax <- str_c("weeks\n", str_c(starts[-(K+1)], ends[-1], sep = "-"))

# define paths
release_folder <- here::here("release20221006")
image_path <- file.path(release_folder, "images")

data_fu <- bind_rows(
  lapply(
    seq_along(subgroups),
    function(x)
    readr::read_csv(file.path(release_folder, glue("check_fu_{x}.csv"))) %>%
      mutate(subgroup=x)
  )
)

plot_fu <- function(.data, periods=12, end_date=study_parameters$end_date) {
  
  plot_data <- .data %>%
    # check that the subgroup x vaccine groups are correct 
    # if not work out if re-extraction is necessary 
    # (might be possible to work it out from the order??)
    group_by(subgroup, k) %>%
    mutate(group = row_number()) %>%
    ungroup() 
  
  if (periods == 12) {
    legend_text_width <- 100
    legend_x_pos <- 0.275
    plot_data <- plot_data %>%
      # remove the BNT162b2 group from subgroup 3 
      filter(!(subgroup == 3 & group == 2))
  } else {
    legend_x_pos <- 0.225
    legend_text_width <- 25
  }
  
  plot_data <- plot_data %>%
    filter(k <= periods, date <= as.Date(end_date)) %>%
    # sum over vaccine groups within subgroup/k/date
    group_by(subgroup, k, date) %>%
    mutate(across(n, sum)) %>%
    ungroup() %>%
    mutate(across(subgroup, factor, labels = subgroups))%>%
    mutate(across(k, factor, levels = 1:K, labels = weeks_since_2nd_vax))
  
  linetype_palette <- c(
    "Delta became dominant" = "solid",
    "First cases of Omicron detected" = "dotted",
    "Omicron became dominant" = "dashed"
  )
  
  names(linetype_palette) <- str_wrap(names(linetype_palette), legend_text_width)
  
  plot_data %>%
    mutate(
      Subgroup=factor(
        subgroup, 
        labels=str_wrap(levels(.$subgroup), width = legend_text_width)
        )
      ) %>%
    ggplot(aes(x = date, y = n, colour = Subgroup, fill = Subgroup, group = glue("{Subgroup}, {k}"))) +
    geom_area(alpha=0.25, position = "identity") +
    # geom_vline(xintercept = c(delta_start, omicron_start, delta_end), linetype = unname(linetype_palette)) +
    geom_vline(aes(xintercept = delta_start, linetype = names(linetype_palette)[1])) +
    geom_vline(aes(xintercept = omicron_start, linetype = names(linetype_palette)[2])) +
    geom_vline(aes(xintercept = delta_end, linetype = names(linetype_palette)[3])) +
    facet_grid(rows = "k") +
    scale_x_date(
      date_breaks = "2 months",
      date_labels = "%b %Y",
      expand = c(0.01,0.01)
    ) +
    scale_linetype_manual(
      name = "Variant dates of interest",
      values = linetype_palette
    ) +
    labs(
      x = "Date of follow-up", 
      y = "Number of individuals followed up (x 1000)"
    ) +
    theme_bw() +
    theme(
      
      axis.title.x = element_text(margin=margin(t=10)),
      axis.title.y = element_text(margin=margin(r=10)),
      
      strip.text.y = element_text(angle = 0),
      strip.background = element_blank(),
      
      legend.position = c(legend_x_pos, 0.2),
      legend.box.background = element_rect(colour = "black", fill = "white")
      
    )
  
}

# plot for the extended follow-up manuscript
plot_fu(.data=data_fu)

ggsave(
  filename = file.path(image_path, glue("fu_all.png")),
  width=15, height=20, units="cm"
)

# combined plot for the 6-month follow-up
plot_fu(.data=data_fu, periods=6, end_date="2021-12-15")

ggsave(
  filename = file.path(image_path, glue("fu_all_6.png")),
  width=15, height=20, units="cm"
)

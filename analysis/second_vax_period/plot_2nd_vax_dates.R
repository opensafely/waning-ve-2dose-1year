################################################################################

# This script:
# - plots and saves the distribution of second vaccination dates

################################################################################

library(tidyverse)
library(glue)
library(lubridate)

################################################################################
# create folder for plots
fs::dir_create(here::here("output", "second_vax_period", "images"))

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
# read data
data_vax_plot <- readr::read_rds(
  here::here("output", "second_vax_period", "data", "data_vax_plot.rds")) %>%
  select(-n) %>%
  # round up to nearest 7 for plot
  mutate(across(n_brand, ~ceiling_any(.x, to=7)))

second_vax_period_dates <- readr::read_rds(
  here::here("output", "second_vax_period", "data", "second_vax_period_dates.rds")) 

# elig_dates info for plot titles
group_age_ranges <- readr::read_csv(
  here::here("output", "lib", "group_age_ranges.csv"))

################################################################################

plot_2nd_vax_dates_fun <- function(
  data, 
  subtitle_string = group_age_ranges
  ) {
  
  jcvi_group <- unique(data$jcvi_group)
  elig_date <- unique(data$elig_date)
  
  # save data for output checking
  capture.output(
    data %>%
      kableExtra::kable("pipe"),
    file = here::here("output",  "second_vax_period", "images", glue("plot_by_region_{jcvi_group}_{elig_date}.txt")),
    append = FALSE
  )
  
  # plot title
  title_string <- glue("Eligible from {elig_date}")
  
  age_range <- data %>%
    distinct(jcvi_group, elig_date) %>%
    left_join(group_age_ranges,
              by = c("jcvi_group", "elig_date")) %>%
    select(age_range) %>% unlist() %>% unname()
  
  jcvi_group_0 <- str_remove(jcvi_group, "^0")
  # age range for plot title
  subtitle_string <- glue("JCVI group {jcvi_group_0}; age range: {age_range} years")
  
  # define breaks for x axis
  x_breaks <- seq(elig_date + weeks(6),
                  elig_date + weeks(20),
                  28) #28 days
  
  # plot histograms by region
  plot_by_region <- ggplot(NULL, aes(x = dose_2)) +
    # overlapping histogram for each brand, binwdith = 1 day
    geom_bar(data = data %>% filter(brand == "ChAdOx1"), 
             aes(y = n_brand, fill = "ChAdOx1"),
             stat = "identity",  alpha = 0.5, width = 1) +
    geom_bar(data = data %>% filter(brand == "BNT162b2"), 
             aes(y = n_brand, fill = "BNT162b2"), 
             stat = "identity", alpha = 0.5, width = 1) +
    # vertical lines to show start and end of second vax period
    geom_vline(data = data, 
               aes(xintercept = start_of_period), colour = "black",
               linetype = "dashed") +
    geom_vline(data = data, 
               aes(xintercept = end_of_period), colour = "black",
               linetype = "dashed") +
    # facet by region
    facet_wrap(~ region, scales = "free_y") +
    scale_x_continuous(breaks = x_breaks,
                       labels = sapply(x_breaks, function(x) str_c(day(x), " ", month(x, label=TRUE)))) +
    scale_y_continuous(expand = expansion(mult = c(0,.05))) +
    scale_fill_discrete(name = "Brand") +
    labs(x = "Date of second vaccination", y = "Number of individuals",
         title = title_string, subtitle = subtitle_string
         # , caption = "Dashed vertical lines show start and end of 28-day second vaccination period."
         ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.box="vertical",
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text.x = element_text(size = 6),
          plot.caption = element_text(size = 10),
          plot.margin = margin(t=0.2, r=0.5, b=0.2, l=0.2, "cm")) 
  # caption:
  # X-axes restricted to 6 to 16 weeks after eligibility date.
  # Bars show the number of individuals who received a second dose of the given brand of vaccine on the given date.
  
  # save the plot
  ggsave(plot_by_region,
         filename = here::here("output", "second_vax_period", "images", glue("plot_by_region_{jcvi_group}_{elig_date}.png")),
         width=20, height=14, units="cm")
  
}

################################################################################

# generate and save plots
lapply(
  data_vax_plot %>%
    left_join(second_vax_period_dates %>%
                select(jcvi_group, elig_date, region, start_of_period, end_of_period), 
              by = c("region", "jcvi_group", "elig_date")) %>% 
    group_split(jcvi_group, elig_date),
       function(x)
         try(plot_2nd_vax_dates_fun(data = x)))


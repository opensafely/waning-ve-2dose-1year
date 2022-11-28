library(tidyverse)

subgroups <- readRDS("~/Documents/waning-ve-2dose-1year/analysis/lib/subgroups.rds")

df_horne <- readr::read_csv("release20221006/estimates_all.csv") %>%
  filter(
    subgroup %in% 1:4,
    outcome %in% c("postest", "covidadmitted"),
    model %in% "adjusted",
    variable %in% "k",
    !reference_row,
    comparison %in% "BNT162b2"
    ) %>%
  rename(k=period) %>%
  mutate(
    period = case_when(
      subgroup %in% 1 ~ "omicron from approx. 35 wks",
      subgroup %in% 2 ~ "omicron from approx. 31 wks",
      subgroup %in% 3 ~ "omicron from approx. 27 wks",
      subgroup %in% 4 ~ "omicron from approx. 15 wks",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(across(subgroup, 
                ~as.character(
                  factor(
                    .x, 
                    labels = subgroups
                    )))) %>%
  mutate(across(outcome, 
                ~as.character(
                  factor(
                    .x, 
                    levels = c("postest", "covidadmitted"),
                    labels = c("infection", "severe")
                    )))) %>%
  mutate(across(c(estimate, conf.low, conf.high),
                ~ 100*(1-exp(.x)))) %>%
  arrange(subgroup, outcome, k) %>%
  mutate(
    fu_lower = (k-1)*28 + 15,
    fu_upper = k*28 + 14,
    fu = str_c(fu_lower, "-", fu_upper)
  ) %>%
  select(age = subgroup, outcome, period,
         ve = estimate, ve_lower = conf.low, ve_upper = conf.high,
         starts_with("fu")) %>%
  mutate(
    author = "horne", 
    country = "England", 
    brand = "BNT", 
    variant = "mix"
  )
  
  

df <- readxl::read_excel(
  "review/waning-beyond-6m.xlsx"
  ) %>%
  filter(
    # must be studying BNT162b2
    str_detect(brand, "BNT"),
    # not only studying HCWs
    !str_detect(age, "HCW")
  )  %>%
  mutate(across(`days fu`, 
                ~ case_when(
                  .x %in% ">180" ~ ">=181",
                  .x %in% "1-3 months" ~ "0-90",
                  .x %in% "4-6 months" ~ "91-180",
                  .x %in% ">=7 months" ~ ">=181",
                  TRUE ~ .x
                ))) %>%
  mutate(period = case_when(
    `end fu` < "2021-12-01" | variant %in% "delta" ~ "non_omciron",
    TRUE ~ "omicron"))


# clean estimates data
ve_split <- str_split(df$`estimated VE`, "\\(")

ve <- sapply(ve_split, function(x) x[1])
ci <- str_split(sapply(ve_split, function(x) x[2]), ",|to")
ve_lower <- sapply(ci, function(x) x[1])
ve_upper <- sapply(ci, function(x) x[2])

f_estimates_clean <- function(x) {
  x <- str_replace_all(x, "·", ".")
  x <- str_replace_all(x, "−", "-")
  x <- str_replace_all(x, "–", "-")
  x <- str_remove_all(x, "\\)")
  x <- str_remove_all(x, "\\s")
  as.numeric(x)
}

estimates_clean <- lapply(
  list(ve, ve_lower, ve_upper),
  f_estimates_clean
)

names(estimates_clean) <- c("ve", "ve_lower", "ve_upper")

# clean fu data

fu_split <- str_split(df$`days fu`, "-")
fu_lower <- as.numeric(str_remove(sapply(fu_split, function(x) x[1]), ">="))
fu_upper <- as.numeric(sapply(fu_split, function(x) x[2]))
fu_clean <- list(fu_lower = fu_lower, fu_upper = fu_upper)


df_clean <- df %>%
  select(author, country, brand, outcome, variant, period, design, age, fu = `days fu`) %>%
  bind_cols(fu_clean, estimates_clean) %>%
  mutate(across(fu_upper, ~if_else(author=="kerr", fu_lower-2.5, .x))) %>%
  mutate(across(fu_upper, ~if_else(author=="kerr", fu_lower+5, .x))) %>%
  mutate(across(age, ~str_c(.x, " years"))) %>%
  bind_rows(df_horne) %>%
  bind_rows(df_horne %>% filter(outcome=="infection") %>% mutate(outcome = "symptomatic infection")) %>%
  mutate(linetype = author=="horne") %>%
  mutate(across(variant, ~if_else(str_detect(.x, "BA"), "omicron", .x))) %>%
  mutate(age_grp = if_else(age %in% c(">=60 years", ">=65 years", ">=70 years", "65+ years"), "older adults", "all ages")) %>%
  mutate(info = str_wrap(str_c(str_to_title(author), ", ", age, ", ", str_to_title(country), ", ", period), 50)) %>%
  mutate(across(fu_upper, ~if_else(is.na(.x), fu_lower  + 30, .x))) %>%
  filter(!(author == "Buchan" & variant == "omicron")) %>%
  mutate(across(outcome, 
                ~ if_else(str_detect(.x, "infection"),
                          str_replace(str_replace(.x, "infection", "SARS CoV-2 infection"), "symptomatic", "Symptomatic"),
                          "Severe COVID-19")
                )) %>%
  filter(author!="kerr") # because only up to 3 months



# outcome_select <- "infection"
# outcome_select <- "symptomatic infection"
# outcome_select <- "severe"
# age_select <- "older adults"
# age_select <- "all ages"

plot_estimates <- function(
  outcome_select, age_select, period_select
) {
  
  # gg plot pallete
  gg_color_hue <- function(n, transparency = 1) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100, alpha = transparency)[1:n]
  }
  
  # lower limit for 7-axis
  lower <- -60
  
  # select data
  plot_data <- df_clean %>%
    filter(
      outcome %in% outcome_select,
      age_grp %in% age_select,
      (period %in% period_select | str_detect(period, "omicron from"))
    ) 
  
  check <- plot_data %>% filter(period %in% period_select)
  if (nrow(check)==0) stop(glue::glue("No {period_select} studies"))
  
  # define palettes
  key <- sort(unique(plot_data$info))
  colour_palette <- gg_color_hue(n = length(key))
  names(colour_palette) <- key
  fill_palette <- gg_color_hue(n = length(key), transparency = 0.2)
  names(fill_palette) <- key
  linetype_palette <- rep("solid", length(key))
  linetype_palette[str_detect(key, "^Horne")] <- "dashed"
  
  plot_data %>%
    mutate(across(starts_with("fu_"), ~.x/7)) %>%
    ggplot() +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_rect(
      aes(
        xmin = fu_lower, xmax = fu_upper, 
        ymin = ve_lower, ymax = ve_upper, 
        fill = info
        ), 
      alpha = 0.2, colour = NA) +
    geom_linerange(
      aes(
        xmin = fu_lower, xmax = fu_upper, 
        y = ve, 
        colour = info, linetype = linetype
        )
      ) +
    scale_y_continuous(limits = c(NA, 100), breaks = seq(lower,100,20)) +
    coord_cartesian(ylim = c(lower,100)) +
    guides(
      fill = guide_legend(
        ncol = 2
      ),
      colour = guide_legend(
        ncol = 2,
        override.aes = list(linetype = linetype_palette, fill = fill_palette)
        )
      ) +
    scale_linetype_discrete(guide = "none") +
    labs(
      title = str_c(outcome_select, ", ", age_select),
      x = "Weeks since second dose",
      y = "Estimated vaccine effectiveness (%)"
      ) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.key.size = unit(1, "cm")
      )
  
  ggsave(
    filename = here::here("review", str_replace_all(glue::glue("{outcome_select}_{age_select}_{period_select}.png"), "\\s", "_")),
    width = 20, height = 14, units = "cm"
  )
  
}

for (o in unique(df_clean$outcome)) {
  for (a in unique(df_clean$age_grp)) {
    for (p in unique(df$period)) {
      try(plot_estimates(
        outcome_select = o,
        age_select = a,
        period_select = p
      ))
    }
  }
}





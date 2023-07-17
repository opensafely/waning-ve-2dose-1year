################################################################################
# This script:
# - combines estimates from all models into csv for release

################################################################################
library(tidyverse)
library(RColorBrewer)
library(lubridate)
library(glue)
library(cowplot)

################################################################################

release_folder <- here::here("release20221006")

fs::dir_create(file.path(release_folder, "images"))

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds")
)
outcomes <- outcomes[outcomes!="covidemergency"]
outcomes <- outcomes[c(3,4,2,5,1)]
outcomes_long <- names(outcomes)
# title case (can't use str_to_title() because it makes caps in words lowercase)
outcomes_long[outcomes=="postest"] <- "Positive Test for SARS-CoV-2"
outcomes_long[outcomes=="covidadmitted"] <- "COVID-19 Hospitalization"
outcomes_long[outcomes=="coviddeath"] <- "COVID-19 Death"
outcomes_long[outcomes=="noncoviddeath"] <- "Non-COVID-19 Death"
names(outcomes) <- outcomes_long
rm(outcomes_long)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_plot_labels <- subgroups
subgroup_plot_labels <- if_else(
  subgroup_plot_labels %in% c("40-64 years", "18-39 years"),
  str_c(subgroup_plot_labels, "*"),
  subgroup_plot_labels
)
subgroup_plot_labels <- str_to_title(subgroup_plot_labels)
subgroup_plot_labels <- str_replace_all(subgroup_plot_labels, "And", "&")
subgroup_plot_labels <- str_replace_all(subgroup_plot_labels, "-", "\u2013")
subgroup_plot_labels <- str_replace_all(subgroup_plot_labels, "65\\+", "\u226565")

subgroup_labels <- seq_along(subgroups)


# define comparisons
comparisons <- c("BNT162b2", "ChAdOx1", "both")

if(Sys.getenv("OPENSAFELY_BACKEND") %in% "") {
  check_fu_path <- release_folder
  surv_path <- release_folder
} else {
  check_fu_path <- "output/tte/images"
  surv_path <- "output/subsequent_vax/tables"
}


################################################################################
# min and max follow-up dates per subgroup

min_max_fu_dates <- bind_rows(lapply(
  1:4,
  function(x)
    readr::read_csv(file.path(check_fu_path, glue("check_fu_{x}.csv"))) %>%
    summarise(across(date, list(min = min, max = max))) %>%
    mutate(subgroup = subgroups[x])
)) %>%
  mutate(across(starts_with("date"), ~format(.x, format = "%d %B %Y")))

caption_str <- min_max_fu_dates %>%
  transmute(text = glue("{date_min} to {date_max} in the \'{subgroup}\' subgroup")) %>%
  unlist() %>% unname()
caption_str <- str_c("Earliest to latest follow-up dates are as follows: ",
                 str_c(caption_str, collapse = "; "),
                 ".",
                 collapse = "")
  
# read estimates data
estimates_all <- readr::read_csv(file.path(release_folder, "estimates_all.csv")) 

# cumulative incidence data
survtable_redacted <- readr::read_csv(
  file.path(surv_path, "survtable_redacted.csv")) %>%
  mutate(across(subgroup,
                ~ case_when(
                  str_detect(.x, "65") ~ 1,
                  str_detect(.x, "18-64") ~ 2,
                  str_detect(.x, "40-64") ~ 3,
                  str_detect(.x, "18-39") ~ 4,
                  TRUE ~ NA_real_
                ))) %>%
  mutate(across(subgroup,
                factor,
                levels = seq_along(subgroups),
                labels = subgroup_plot_labels
  )) 

################################################################################
# gg plot pallete
# gg_color_hue <- function(n, transparency = 1) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100, alpha = transparency)[1:n]
# }

formatpercent100 <- function(x,accuracy){
  formatx <- scales::label_percent(accuracy)(x)
  
  if_else(
    formatx==scales::label_percent(accuracy)(1),
    paste0(">",scales::label_percent(1)((100-accuracy)/100)),
    formatx
  )
}

# scale for x-axis
K <- study_parameters$K
ends <- seq(2, (K+1)*4, 4)
starts <- ends + 1
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "\u2013")
weeks_since_2nd_vax_odd <- weeks_since_2nd_vax
weeks_since_2nd_vax_odd[seq(2, length(weeks_since_2nd_vax), 2)] <- ""

################################################################################

# y_lab <- "Hazard Ratio (HR)"
# y_lab_2 <- "Estimated vaccine effectiveness = 100*(1-HR)\n "
# y_lab_adj <- "adjusted Hazard Ratio (aHR)"
# y_lab_adj_2 <- "Estimated vaccine effectiveness = 100*(1-aHR)\n "
# x_lab <- "Weeks since second dose"
# # x_lab_nofootnote <- str_remove(x_lab, "\\*\\*")

legend_width <- 15

plot_data <- estimates_all %>%
  filter(
    !reference_row,
    variable %in% "k"
  ) %>%
  # keep only outcomes of interest
  filter(outcome %in% outcomes) %>%
  mutate(across(c(estimate, conf.low, conf.high), exp)) %>%
  mutate(k=as.integer(label)) %>%
  mutate(across(subgroup,
                factor,
                levels = subgroup_labels,
                labels = subgroups
  )) %>%
  mutate(across(model,
                factor,
                levels = c("unadjusted", "adjusted"),
                labels = sapply(c("Stratfied Cox model, no further adjustment", 
                                  "Stratfied Cox model, adjustment for demographic and clinical variables"),
                                str_wrap, width=100))) %>%
  mutate(outcome_unlabelled = outcome) %>%
  mutate(across(outcome,
                factor,
                levels = unname(outcomes),
                labels = str_wrap(names(outcomes), 17)
  )) %>%
  mutate(across(subgroup,
                factor,
                levels = subgroups,
                labels = subgroup_plot_labels
  )) %>%
  mutate(line_group = str_c(subgroup, comparison, outcome,  model, sep = "; "))


# spacing of points on plot
position_dodge_val <- 0.6

# colours of points
palette_adj <- RColorBrewer::brewer.pal(3, "Dark2")
names(palette_adj) <- c(comparisons[1:2], "Unvaccinated")

# # shapes of points
variants <- c("alpha", "delta", "omicron")
variant_shapes <- c(21,24,22)
names(variant_shapes) <- variants

# fill of points
fill_groups <- c("other_BNT162b2", "other_ChAdOx1", "omicron_BNT162b2", "omicron_ChAdOx1")
palette_fill <- c("white", "white", palette_adj[1:2])
names(palette_fill) <- fill_groups

# breaks and lims for y-axes
primary_vax_y1 <- lst(breaks = c(0.005, 0.05, 0.25, 1, 4), 
                       limits = range(breaks))

# cumulative incidence (plot A)

# scale for x-axis
x_breaks <- seq(3, 48, 4)
x_labels <- as.character(x_breaks)
alpha_area <- 0.5
# create plot
# max_nrisk <- 1000000#max(survtable_redacted$n.risk)

plot_ci_data <- survtable_redacted %>%
  filter(
    !(as.integer(subgroup) == 3L & arm == "BNT162b2"),
    !(as.integer(subgroup) == 4L & arm == "ChAdOx1")
  ) %>%
  filter(time <= 48) %>%
  # rescale time to "weeks since second dose" - it's currently weeks since start of comparison period 1
  mutate(time = time + 2) %>%
  mutate(
    variant = factor(
      case_when(
        as.integer(subgroup) == 1L & time <= 1*4 ~ "alpha",
        as.integer(subgroup) == 1L & time < 9*4 ~ "delta",
        as.integer(subgroup) == 1L ~ "omicron",
        as.integer(subgroup) == 2L & time < 8*4 ~ "delta",
        as.integer(subgroup) == 2L ~ "omicron",
        as.integer(subgroup) == 3L & time < 7*4 ~ "delta",
        as.integer(subgroup) == 3L ~ "omicron",
        as.integer(subgroup) == 4L & time < 4*4 ~ "delta",
        as.integer(subgroup) == 4L ~ "omicron",
        TRUE ~ NA_character_
      ),
      levels = variants
    )
  ) %>%
  mutate(
    c.inc_alphadelta = if_else(variant %in% c("alpha", "delta"), c.inc, NA_real_),
    c.inc_omicron = if_else(variant %in% c("omicron"), c.inc, NA_real_)
  ) 

# get earliest time when ci over 80%
ci_80_time <- plot_ci_data %>%
  arrange(subgroup, arm, time) %>%
  filter(c.inc>0.8) %>%
  group_by(subgroup, arm) %>%
  slice(1) %>%
  select(subgroup, arm, time) %>%
  ungroup() %>%
  mutate(
    time_min = round(time,0),
    # plus 1 because people followed up to end of the period before excluded for subsequent dose
    k_min = as.integer(cut(time_min, breaks = seq(2, (K+1)*4, 4))) + 1
    ) 
  # use distinct to avoid overlapping rectangles
  # this is ok because the value is the same across arms
  # distinct(subgroup, k_min, .keep_all = TRUE)

plot_ci <- plot_ci_data %>%
  # add time_min for shaded rectangles
  bind_rows(ci_80_time) %>%
  ggplot() +
  geom_rect(
    aes(
      xmin = time_min, 
      xmax = 49,
      ymin = 0,
      ymax = 1
    ), 
    fill = "grey70",
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = time,
      colour = arm,
      y = c.inc_alphadelta
      ),
    linetype = "dashed",
  ) +
  geom_line(
    aes(
      x = time,
      colour = arm,
      y = c.inc_omicron
      ),
  ) +
  facet_grid(
    . ~ subgroup,
    switch = "y",
    space = "free_x"
  ) +
  scale_x_continuous(
    name = "Weeks Since Second Dose",
    expand = expansion(add = c(1, 0)),
    limits = c(3,49), # set max as 49 to align with panel B, which is plotted at midpoint in [47,50] interval, which is 49
    breaks = seq(3,50,4)#,
    # labels = NULL
  ) +
  scale_y_continuous(
    name = str_wrap("Cumulative Incidence of Subsequent Dose", 14),
    expand = c(0,0),
    limits = c(0,1),
    labels = format(seq(0,1,0.25), nsmall=2),
    oob = scales::oob_keep
  ) +
  scale_colour_manual(
    name = NULL,
    values = palette_adj,
    guide = "none"
  ) +
  labs(x = NULL) +
  theme_bw(base_family = "helvetica") +
  theme(
    
    axis.line = element_blank(),
    # axis.line = element_line(colour = "black", linewidth = 0.35), # 0.35mm ~= 1pt
    
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(
      size=8,
      angle = 45
    ),
    
    axis.title.x = element_text(
      size = 8, 
      margin = margin(t = 5, r = 0, b = 5, l = 0)
    ),
    axis.title.y.left = element_text(
      size = 8, 
      margin = margin(t = 0, r = 10, b = 0, l = 0),
      angle = 90,
      vjust = 0.5
    ),
    axis.title.y.right = element_text(
      size = 8,
      margin = margin(t = 0, r = 0, b = 0, l = 10),
      angle = 0,
      vjust = 0.5
    ),
    
    axis.ticks = element_line(linewidth = 0.35),
    
    panel.border = element_rect(linewidth = 0.35, colour = "black"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size=8),
    
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    plot.margin = margin(l = 20, t = 5, r = 5, b = 5),
    
  )

# vaccine vs unvaccinated (plot B)
plot_vax_data <- plot_data %>%
  filter(
    comparison != "both",
    outcome_unlabelled != "anytest",
    as.integer(model) == 2,
    !(as.integer(subgroup) == 3L & comparison == "BNT162b2"),
    !(as.integer(subgroup) == 4L & comparison == "ChAdOx1")
  ) %>%
  mutate(
    across(
      outcome,
      factor,
      levels = levels(plot_data$outcome)#,
      # labels = str_c("HR for\n", levels(plot_data$outcome))
      )
    ) %>%
  droplevels() %>%
  mutate(
    variant = factor(
      case_when(
        as.integer(subgroup) == 1L & period == 1 ~ "alpha",
        as.integer(subgroup) == 1L & period < 9 ~ "delta",
        as.integer(subgroup) == 1L ~ "omicron",
        as.integer(subgroup) == 2L & period < 8 ~ "delta",
        as.integer(subgroup) == 2L ~ "omicron",
        as.integer(subgroup) == 3L & period < 7 ~ "delta",
        as.integer(subgroup) == 3L ~ "omicron",
        as.integer(subgroup) == 4L & period < 4 ~ "delta",
        as.integer(subgroup) == 4L ~ "omicron",
        TRUE ~ NA_character_
      ),
      levels = variants
    )
  ) %>%
  mutate(
    fill_group = factor(
      str_c(
        if_else(variant %in% "omicron", "omicron", "other"),
        comparison,
        sep = "_"
      ),
      levels = fill_groups
    )
  ) 


plot_vax <- plot_vax_data %>%
  bind_rows(
    plot_vax_data %>% 
      distinct(subgroup) %>%
      mutate(subgroup_join=as.integer(subgroup)) %>%
      inner_join(
        ci_80_time %>% 
          mutate(subgroup_join = as.integer(subgroup)) %>%
          mutate(
            comparison = factor(
              arm, levels = levels(plot_vax_data$comparison)
              )
            ) %>%
          select(subgroup_join, comparison, k_min),
        by = join_by(subgroup_join)
      ) %>%
      select(-subgroup_join) %>%
      expand_grid(
        outcome = levels(plot_vax_data$outcome)
        )
  ) %>%
  # add a dummy row so that the unvaccinated line is shown on the legend
  add_row(
    comparison = "Unvaccinated",
    outcome = factor(levels(plot_vax_data$outcome)[1], levels = levels(plot_vax_data$outcome)),
    subgroup = factor(levels(plot_vax_data$subgroup)[1], levels = levels(plot_vax_data$subgroup))
    ) %>%
  mutate(across(comparison, factor, levels = names(palette_adj))) %>%
  ggplot() +
  # add shaded areas where cumulative incidence of third dose over 80%
  geom_rect(
    aes(
      xmin = k_min - 0.5, 
      xmax = 12 + 0.5,
      ymin = primary_vax_y1$limits[1],
      ymax = primary_vax_y1$limits[2]
      ), 
    fill = "grey70",
    alpha = 0.2
  ) +
  geom_hline(
    aes(yintercept=1),
    colour='black', linetype = "dashed", linewidth = 0.35
    ) +
  geom_linerange(
    aes(
      ymin = conf.low, ymax = conf.high,
      x = k, 
      colour = comparison, 
      shape = variant,
      # alpha = variant,
      fill = fill_group
      ),
    position = position_dodge(width = position_dodge_val)
  ) +
  geom_point(
    aes(
      y = estimate,
      x = k, 
      colour = comparison, 
      shape = variant,
      # alpha = variant,
      fill = fill_group
      ),
    position = position_dodge(width = position_dodge_val),
    size = 0.9,
    shape = 21
  ) +
  # lemon::facet_rep_grid(
    # repeat.tick.labels = TRUE,
  facet_grid(
    outcome ~ subgroup, 
    switch = "y",
    # scales = "free_x",
    space = "free_x"
    ) +
  labs(x = "Weeks Since Second Dose (4-week Comparison Period)") +
  scale_x_continuous(
    expand = c(0,0),
    breaks = 1:12, # scale is time since start of period 1
    labels = weeks_since_2nd_vax # label scale as time since second vax
  ) +
  scale_y_log10(
    expand = c(0,0),
    # expand = expansion(mult = c(0,0.05)),
    name = "Hazard Ratio",
    breaks = primary_vax_y1[["breaks"]],
    limits = primary_vax_y1[["limits"]],
    oob = scales::oob_keep
  ) +
  scale_fill_discrete(guide = "none") +
  scale_color_manual(
    values = palette_adj,
    labels = names(palette_adj),
    drop = FALSE,
    name = NULL,
    guide = "none"
    ) +
  scale_fill_manual(
    values = palette_fill, 
    name = NULL,
    guide = "none"
  ) +
  scale_linetype_manual(
    values = c("dashed", "solid"), 
    name = NULL, 
    guide = "none"
    ) +
  guides(
    colour = guide_legend(
      title = expression(underline(Vaccination~Group)),
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        shape = NA,
        size = 3
      )
    )
    ) +
  theme_bw(base_family = "helvetica") +
  theme(
    
    axis.line = element_blank(),
    # axis.line = element_line(colour = "black", linewidth = 0.35), # 0.35mm ~= 1pt
    
    axis.title.x = element_text(
      size=8, 
      margin = margin(t = 5, r = 0, b = 0, l = 0)
      ),
    axis.title.y = element_text(
      size = 8
      ),
    
    axis.text.x = element_text(size=8, angle=45, hjust=1),
    axis.text.y = element_text(size=8),
    axis.ticks = element_line(linewidth = 0.35),
    
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_blank(),
    strip.text.y.left = element_text(angle = 90),
    strip.text = element_text(size=8),
    
    panel.border = element_rect(linewidth = 0.35, colour = "black"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    plot.margin = margin(l = 20, t = 10, r = 5, b = 5),
    
    legend.position = c(0.875, 0.63),
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size=8),
    legend.text = element_text(size=8)
    # legend.background = element_rect(linetype = 1, linewidth = 0.35, colour = "black"),
    
  ) 

set_null_device(cairo_pdf)
plot_combined <- plot_grid(
  plot_ci, plot_vax, 
  nrow = 2, rel_heights = c(0.23,0.77),
  labels = c("A)", "B)"),
  label_fontface = "plain",
  label_fontfamily = "helvetica",
  label_size = 8,
  align="v", axis = c("lr")
)

print(caption_str)

# save the plot
ggsave(plot_combined,
       filename = file.path(release_folder, "images", glue("hr_vax_ci.pdf")),
       device = cairo_pdf,
       width=9, height=6, units="in"
       )

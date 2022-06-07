################################################################################
# This script:
# - combines estimates from all models into csv for release

################################################################################
library(tidyverse)
library(RColorBrewer)
library(lubridate)
library(glue)

################################################################################
if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds")
)
outcomes <- outcomes[outcomes!="covidemergency"]
outcomes_order <- c(3,4,2,5,1)
outcomes_long <- names(outcomes)
outcomes_long[outcomes=="covidadmitted"] <- "COVID-19 hospitalisation"
names(outcomes) <- outcomes_long
rm(outcomes_long)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# define comparisons
comparisons <- c("BNT162b2", "ChAdOx1", "both")

################################################################################
# min and max follow-up dates per subgroup
min_max_fu_dates <- readr::read_csv(
  here::here(release_folder, glue("data_min_max_fu.csv"))) %>%
  mutate(across(ends_with("date"),
                ~ str_c(day(.x), " ", month(.x, label=TRUE))))

# read estimates data
estimates_all <- 
  bind_rows(
    readr::read_csv(here::here(release_folder, "estimates_all.csv")),
    readr::read_csv(here::here(release_folder, "estimates_6575.csv"))
  ) %>%
  mutate(
    sex = if_else(
      str_detect(subgroup, "Female|Male"),
      str_extract(subgroup, "Female|Male"),
      "Both"
    ),
    ageband = case_when(
      str_detect(subgroup, "_65") ~ "65-74 years",
      str_detect(subgroup, "_75") ~ "75+ years",
      TRUE ~ "all"
    ),
    subgroup = as.integer(str_extract(subgroup, "^\\d"))
  )

# read metareg data
metareg_results_k <- readr::read_rds(
  here::here(release_folder, "metareg_results_k.rds")) %>%
  select(subgroup, comparison, sex, ageband, outcome, k, starts_with("line")) %>%
  mutate(model="adjusted") 

################################################################################
# gg plot pallete
gg_color_hue <- function(n, transparency = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = transparency)[1:n]
}

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
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")

################################################################################

page_height <- 27
page_width <- 16

y_lab <- "Hazard Ratio (HR)"
y_lab_2 <- "Estimated vaccine effectiveness = 100*(1-HR)\n "
y_lab_adj <- "adjusted Hazard Ratio (aHR)"
y_lab_adj_2 <- "Estimated vaccine effectiveness = 100*(1-aHR)\n "
x_lab <- "Weeks since second dose**"
# x_lab_nofootnote <- str_remove(x_lab, "\\*\\*")

legend_width <- 15
subgroup_plot_labels <- if_else(
  subgroups %in% c("40-64 years", "18-39 years"),
  str_c(subgroups, "*"),
  subgroups
)

plot_data <- estimates_all %>%
  filter(
    !reference_row,
    variable %in% "k"
  ) %>%
  # keep only outcomes of interest
  filter(outcome %in% outcomes) %>%
  # remove as too few events
  filter(
    !(comparison %in% c("BNT162b2", "both") & subgroup %in% c(3,4) & outcome == "noncoviddeath"),
    !(comparison %in% c("BNT162b2", "both") & subgroup %in% 3 & outcome == "covidadmitted")
  ) %>%
  mutate(across(c(estimate, conf.low, conf.high), exp)) %>%
  mutate(k=as.integer(label)) %>%
  left_join(
    min_max_fu_dates %>%
      left_join(
        tibble(subgroup = subgroups, subgroup_label = subgroup_labels),
        by = "subgroup"
      ) %>%
      select(-subgroup) %>%
      rename(subgroup = subgroup_label), 
    by = "subgroup"
  ) %>%
  mutate(order1 = 10*k) %>%
  mutate(k_labelled = k) %>%
  mutate(across(k_labelled, 
                factor, 
                levels = 1:K,
                labels = weeks_since_2nd_vax)) %>%
  mutate(across(subgroup,
                factor,
                levels = subgroup_labels,
                labels = subgroups
  )) %>%
  mutate(k_labelled_dates = k_labelled) %>%
  mutate(across(k_labelled_dates,
                ~ case_when(
                  k == 1
                  ~ str_c(.x, "\n(from\n", min_fu_date, ")"),
                  k == K
                  ~ str_c(.x, "\n(to\n", max_fu_date, ")"),
                  TRUE ~ as.character(.x)))) %>%
  arrange(k) %>%
  group_by(k, subgroup, sex, ageband, outcome, comparison) %>%
  mutate(order2 = row_number()) %>%
  ungroup() %>%
  mutate(order = order1 + order2) %>%
  left_join(
    metareg_results_k, 
    by = c("subgroup", "comparison", "sex", "ageband", "outcome", "model", "k")
  ) %>%
  mutate(across(model,
                factor,
                levels = c("unadjusted", "adjusted"),
                labels = sapply(c("Stratfied Cox model, no further adjustment", 
                                  "Stratfied Cox model, adjustment for demographic and clinical variables"),
                                str_wrap, width=100))) %>%
  mutate(outcome_unlabelled = outcome) %>%
  mutate(across(outcome,
                factor,
                levels = unname(outcomes[outcomes_order]),
                labels = str_wrap(names(outcomes[outcomes_order]), 10)
  )) %>%
  mutate(across(subgroup,
                factor,
                levels = subgroups,
                labels = str_wrap(subgroup_plot_labels, 25)
  )) %>%
  mutate(line_group = str_c(subgroup, sex, ageband, comparison, outcome,  model, sep = "; ")) %>%
  # only plot line within range of estimates
  mutate(k_nonmiss = if_else(!is.na(estimate), k, NA_integer_)) %>%
  group_by(line_group) %>%
  mutate(keep = 0 < sum(!is.na(k_nonmiss))) %>%
  filter(keep) %>%
  mutate(
    min_k = min(k_nonmiss, na.rm = TRUE),
    max_k = max(k_nonmiss, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(across(line,
                ~ if_else(min_k <= k & k <= max_k,
                          .x,
                          NA_real_)))


# spacing of points on plot
position_dodge_val <- 0.6

# shape of points
comparison_shapes <- c(16,17,15)
names(comparison_shapes) <- comparisons
hollow_shapes <- c(21,24,22)
names(hollow_shapes) <- comparisons

comparison_linetypes <- c("dashed", "dotted", "dotdash")
names(comparison_linetypes) <- comparisons

# colours of points
palette_adj <- c(gg_color_hue(n=2), gg_color_hue(n=3)[2])
names(palette_adj) <- comparisons
palette_unadj <- c(gg_color_hue(n=2, transparency=0.3),
                   gg_color_hue(n=3, transparency=0.3)[2])
names(palette_adj) <- comparisons


# point_shapes <- c(22,24)
# breaks and lims for y-axes
primary_vax_y1 <- list(breaks = c(0.02, 0.05, 0.2, 0.5, 1, 2), 
                       limits = c(0.02, 2))
primary_vax_y2 <- list(breaks = c(0,0.5,0.8, 0.95, 0.98))
primary_brand_y1 <- list(breaks = c(0.2, 0.5, 1, 2, 5), 
                         limits = c(0.2, 5))
anytest_y1 <- list(breaks = c(0.5, 1, 2, 5), 
                   limits = c(0.5, 5))

################################################################################
# vaccine vs unvaccinated
plot_vax <- plot_data %>%
  filter(
    sex == "Both",
    comparison != "both",
    ageband == "all",
    outcome_unlabelled != "anytest",
    as.integer(model) == 2
  ) %>%
  droplevels() %>%
  ggplot(aes(
    x = reorder(k_labelled_dates, order), 
    colour = comparison, 
    shape = comparison,
    fill = comparison
  )) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_line(
    aes(y = line, 
        colour = comparison, 
        linetype = comparison,
        group = line_group), 
    alpha = 0.6
  ) +
  geom_linerange(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = position_dodge_val)
  ) +
  geom_point(
    aes(y = estimate),
    position = position_dodge(width = position_dodge_val)
  ) +
  facet_grid(outcome ~ subgroup, switch = "y", scales = "free", space = "free_x") +
  scale_y_log10(
    name = y_lab_adj,
    breaks = primary_vax_y1[["breaks"]],
    limits = primary_vax_y1[["limits"]],
    oob = scales::oob_keep,
    sec.axis = sec_axis(
      ~(1-.),
      name=y_lab_adj_2,
      breaks = primary_vax_y2[["breaks"]],
      labels = function(x){formatpercent100(x, 1)}
    )
  ) +
  labs(
    x = x_lab
  ) +
  scale_fill_discrete(guide = "none") +
  scale_shape_manual(values = comparison_shapes[1:2], name = NULL) +
  scale_color_manual(values = palette_adj[1:2], name = NULL) +
  scale_linetype_manual(values = comparison_linetypes[1:2], name = NULL) +
  guides(shape = guide_legend(
    title = NULL, 
    override.aes = list(colour = palette_adj[1:2], fill = comparison_shapes[1:2])
  )) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=10),
    
    axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size=8),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = c(0.88, 0.14),
    # big margins to cover up grid lines
    legend.margin = margin(t = 30, r = 20, b = 30, l = 10),
    legend.key.width = unit(2, 'cm'),
    # legend.position = "bottom",
    legend.text = element_text(size=10)
  ) 

# save the plot
ggsave(plot_vax,
       filename = here::here(release_folder, glue("hr_vax.png")),
       width=page_height, height=page_width, units="cm")


ggsave(plot_vax + theme(plot.margin = margin(2, 2, 2, 2, "cm")),
       filename = here::here(release_folder, glue("hr_vax.pdf")),
       width=page_height, height=page_width, units="cm")

################################################################################
# brand comparison
plot_brand <- plot_data %>%
  filter(
    sex == "Both",
    ageband == "all",
    comparison == "both",
    outcome_unlabelled != "anytest",
    as.integer(model) == 2
  ) %>%
  droplevels() %>%
  # complete(subgroup, comparison, outcome, k_labelled) %>%
  ggplot(aes(
    x = reorder(k_labelled_dates, order)
  )) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_line(
    aes(y = line, 
        group = line_group), 
    linetype = comparison_linetypes[3],
    colour=palette_adj[3],
    alpha = 0.6
  ) +
  geom_linerange(
    aes(ymin = conf.low, ymax = conf.high), 
    position = position_dodge(width = position_dodge_val),
    color = palette_adj[3],
    fill = palette_adj[3]
  ) +
  geom_point(
    aes(y = estimate),
    position = position_dodge(width = position_dodge_val),
    color = palette_adj[3],
    fill = palette_adj[3],
    shape = comparison_shapes[3]
  ) +
  facet_grid(outcome ~ subgroup, switch = "y", scales = "free", space = "free_x") +
  scale_y_log10(
    name = y_lab_adj,
    breaks = primary_brand_y1[["breaks"]],
    limits = primary_brand_y1[["limits"]],
    oob = scales::oob_keep) +
  guides(
    color = guide_legend(override.aes = list(linetype = 0))
  ) +
  labs(
    x = x_lab
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=10),
    
    axis.title.x = element_text(size = 10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size=8),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic")
    
  ) 
ggsave(plot_brand,
       filename = here::here(release_folder, glue("hr_brand.png")),
       width=page_height, height=page_width, units="cm")
ggsave(plot_brand + theme(plot.margin = margin(2, 2, 2, 2, "cm")),
       filename = here::here(release_folder, glue("hr_brand.pdf")),
       width=page_height, height=page_width, units="cm")

################################################################################
# anytest
# vaccine vs unvaccinated
plot_vax_anytest <- plot_data %>%
  filter(
    sex == "Both",
    ageband == "all",
    comparison != "both",
    outcome_unlabelled == "anytest",
    as.integer(model) == 2
  ) %>%
  ggplot(aes(
    x = reorder(k_labelled_dates, order), 
    y = estimate, 
    colour = comparison, 
    shape = comparison,
    fill = comparison)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_line(
    aes(y = line, 
        colour=comparison, 
        linetype = comparison,
        group = line_group), 
    alpha = 0.6
  ) +
  geom_linerange(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = position_dodge_val)) +
  geom_point(
    position = position_dodge(width = position_dodge_val)
  ) +
  facet_wrap( ~ subgroup, scales = "free") +
  scale_y_log10(
    name = y_lab_adj,
    breaks = anytest_y1[["breaks"]],
    limits = anytest_y1[["limits"]],
    oob = scales::oob_keep
  ) +
  labs(
    x = x_lab,
    title = "Any SARS-CoV-2 test"
  ) +
  scale_shape_manual(values = comparison_shapes[1:2], name = NULL) +
  scale_color_manual(values = palette_adj[1:2], name = NULL) +
  scale_fill_discrete(guide = "none") +
  scale_linetype_manual(values = comparison_linetypes[1:2], name = NULL) +
  guides(shape = guide_legend(
    title = NULL, 
    override.aes = list(colour = palette_adj[1:2], fill = comparison_shapes[1:2])
  )) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=10),
    
    axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size=8),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0, size=10),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = "bottom",
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size=10)
  ) 

# save the plot
ggsave(plot_vax_anytest,
       filename = here::here(release_folder, glue("hr_vax_anytest.png")),
       width=page_width, height=12, units="cm")

################################################################################
# anytest
# brand comparison
plot_brand_anytest <- plot_data %>%
  filter(
    sex == "Both",
    ageband == "all",
    comparison == "both",
    outcome_unlabelled == "anytest",
    as.integer(model) == 2
  ) %>%
  droplevels() %>%
  ggplot(aes(
    x = reorder(k_labelled, order), 
    y = estimate)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_line(
    aes(y = line, 
        group = line_group), 
    colour = palette_adj[3],
    linetype=comparison_linetypes[3], 
    alpha = 0.6
  ) +
  geom_linerange(
    aes(ymin = conf.low, ymax = conf.high), 
    position = position_dodge(width = position_dodge_val),
    color = palette_adj[3],
    fill = palette_adj[3]
  ) +
  geom_point(
    position = position_dodge(width = position_dodge_val),
    color = palette_adj[3],
    fill = palette_adj[3],
    shape = comparison_shapes[3]
  ) +
  facet_wrap( ~ subgroup, scales = "free", ncol=2) +
  scale_y_log10(
    name = y_lab_adj,
    breaks = anytest_y1[["breaks"]],
    limits = anytest_y1[["limits"]],
    oob = scales::oob_keep
  ) +
  labs(
    x = x_lab,
    title = "Any SARS-CoV-2 test"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=10),
    
    axis.title.x = element_text(size = 10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size=8),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0, size=10),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic")
  ) 
ggsave(plot_brand_anytest,
       filename = here::here(release_folder, glue("hr_brand_anytest.png")),
       width=page_width, height=12, units="cm")

################################################################################
# unadjusted and adjusted estimates for each comparison

plot_unadj_adj <- function(plot_comparison) {
  
  i <- which(comparisons == plot_comparison)
  palette <- unname(c(palette_unadj[i], palette_adj[i]))
  names(palette) <- levels(plot_data$model)
  point_shapes <- unname(comparison_shapes[i])
  
  # vaccine vs unvaccinated
  plot_vax_0 <- plot_data %>%
    filter(
      sex == "Both",
      ageband == "all",
      comparison == plot_comparison,
      outcome_unlabelled != "anytest"
    ) %>%
    ggplot(aes(
      x = reorder(k_labelled_dates, order), 
      y = estimate, 
      colour = model, 
      fill = model)) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_linerange(
      aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = position_dodge_val)) +
    geom_point(
      position = position_dodge(width = position_dodge_val),
      shape = point_shapes
    ) +
    facet_grid(outcome ~ subgroup, switch = "y", scales = "free", space = "free_x") 
  
  if (plot_comparison == "both") {
    
    plot_vax_1 <- plot_vax_0  +
      scale_y_log10(
        name = "Hazard ratio (HR)",
        breaks = c(0.02, 0.05, 0.2, 0.5, 1, 2),
        limits = c(0.02, 2),
        oob = scales::oob_keep
      )
    
  } else {
    
    plot_vax_1 <- plot_vax_0  +
      scale_y_log10(
        name = y_lab,
        breaks = primary_vax_y1[["breaks"]],
        limits = primary_vax_y1[["limits"]],
        oob = scales::oob_keep,
        sec.axis = sec_axis(
          ~(1-.),
          name=y_lab_2,
          breaks = c(0,0.5,0.8, 0.95, 0.98),
          labels = function(x){formatpercent100(x, 1)}
        )
      )
    
  }
  
  plot_vax_2 <- plot_vax_1 +
    labs(
      x = x_lab
    ) +
    scale_colour_manual(name = NULL,
                        values = palette) +
    scale_fill_manual(guide = "none",
                      values = palette) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=10),
      
      axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      strip.text = element_text(size=8),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = "bottom",
      legend.text = element_text(size=8)
    ) 
  
  # save the plot
  ggsave(plot_vax_2,
         filename = here::here(release_folder, glue("hr_vax_{plot_comparison}_unadj.png")),
         width=page_height, height=page_width, units="cm")
  
  # plot comparison
  plot_vax_anytest <- plot_data %>%
    filter(
      sex == "Both",
      ageband == "all",
      comparison == plot_comparison,
      outcome_unlabelled == "anytest"
    ) %>%
    ggplot(aes(
      x = reorder(k_labelled_dates, order), 
      y = estimate, 
      colour = model, 
      fill = model)) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_linerange(
      aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = position_dodge_val)) +
    geom_point(
      position = position_dodge(width = position_dodge_val),
      shape = point_shapes
    ) +
    facet_wrap( ~ subgroup, scales = "free", nrow=2) +
    scale_y_log10(
      name = y_lab,
      breaks = anytest_y1[["breaks"]],
      limits = anytest_y1[["limits"]],
      oob = scales::oob_keep
    ) +
    labs(
      x = x_lab,
      title = "Any SARS-CoV-2 test"
    ) +
    scale_colour_manual(name = NULL,
                        values = palette) +
    scale_fill_manual(guide = "none",
                      values = palette) +
    guides(colour=guide_legend(nrow=2)) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=10),
      
      axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      strip.text = element_text(size=8),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0, size=10),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = "bottom",
      legend.text = element_text(size=10)
    ) 
  
  # save the plot
  ggsave(plot_vax_anytest,
         filename = here::here(release_folder, glue("hr_vax_anytest_{plot_comparison}_unadj.png")),
         width=page_width, height=12, units="cm")
  
}

for (i in c("BNT162b2", "ChAdOx1", "both")) {
  plot_unadj_adj(i)
}

################################################################################
# sex stratified results for each comparison

plot_strata <- function(plot_comparison, strata) {
  
  plot_data_0 <- plot_data %>%
    filter(
      comparison == plot_comparison,
      as.integer(model) == 2
    )
  
  leg_pos <- "bottom"
  if (strata == "sex") {
    strata_labs <- c("Female", "Male")
    plot_data_1 <- plot_data_0 %>%
      filter(
        sex != "Both"
      ) %>%
      rename(strata_var = !! strata)
    leg_rows <- 1
    strata_page_width <- 27
    strata_page_height <- 16
  } else if( strata == "ageband") {
    strata_labs <- c("65-74 years", "75+ years")
    plot_data_1 <- plot_data_0 %>%
      filter(
        ageband != "all"
      ) %>%
      rename(strata_var = !! strata)
    leg_rows <- 2
    strata_page_height <- 16
    strata_page_width <- 14
  } else {
    stop("strata must be sex or ageband")
  }
  
  i <- which(comparisons == plot_comparison)
  palette <- palette_adj[i]
  point_shapes <- unname(c(comparison_shapes[i], hollow_shapes[i]))
  names(point_shapes) <- strata_labs
  fill_shapes <- c(palette, "white")
  names(fill_shapes) <- strata_labs
  line_types <- c("dashed", "dotted")
  names(line_types) <- strata_labs
  
  # vaccine vs unvaccinated
  plot_vax_0 <- plot_data_1 %>%
    filter(
      outcome_unlabelled != "anytest"
    ) %>%
    ggplot(aes(
      x = reorder(k_labelled_dates, order), 
      y = estimate,
      colour = comparison,
      fill = strata_var,
      shape = strata_var)) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_line(
      aes(y = line, 
          linetype = strata_var,
          group = line_group), 
      colour = palette,
      alpha = 0.6
    ) +
    geom_linerange(
      aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = position_dodge_val)) +
    geom_point(
      position = position_dodge(width = position_dodge_val)
    ) +
    facet_grid(outcome ~ subgroup, switch = "y", scales = "free", space = "free_x") 
  
  if (plot_comparison == "both") {
    
    plot_vax_1 <- plot_vax_0  +
      scale_y_log10(
        name = y_lab_adj,
        breaks = c(0.02, 0.05, 0.2, 0.5, 1, 2),
        limits = c(0.02, 2),
        oob = scales::oob_keep
      )
    
  } else {
    
    plot_vax_1 <- plot_vax_0  +
      scale_y_log10(
        name = y_lab_adj,
        breaks = primary_vax_y1[["breaks"]],
        limits = primary_vax_y1[["limits"]],
        oob = scales::oob_keep,
        sec.axis = sec_axis(
          ~(1-.),
          name=y_lab_adj_2,
          breaks = c(0,0.5,0.8, 0.95, 0.98),
          labels = function(x){formatpercent100(x, 1)}
        )
      )
    
  }
  
  plot_vax_2 <- plot_vax_1 +
    labs(
      x = x_lab
    ) +
    scale_linetype_manual(
      name = NULL,
      values = line_types
    ) +
    scale_shape_manual(
      values = point_shapes
    ) +
    scale_colour_manual(guide = "none",
                        values = palette) +
    scale_fill_manual(guide = "none",
                      values = fill_shapes) +
    guides(
      shape = guide_legend(
        title = NULL, 
        override.aes = list(colour = palette, fill = fill_shapes),
        nrow=leg_rows, byrow=TRUE
      )
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=10),
      
      axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      strip.text = element_text(size=8),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = leg_pos,
      legend.key.width = unit(2, 'cm'),
      legend.text = element_text(size=8)
    ) 
  
  # save the plot
  ggsave(plot_vax_2,
         filename = here::here(release_folder, glue("hr_vax_{plot_comparison}_{strata}.png")),
         width=strata_page_width, height=strata_page_height, units="cm")
  
  
  strata_page_width <- 16
  strata_page_height <- 12
  if (strata=="ageband") {
    strata_page_height <- 8
    leg_pos <- "right"
  }
  # plot comparison anytest
  plot_vax_anytest <- plot_data_1 %>%
    filter(
      outcome_unlabelled == "anytest"
    ) %>%
    ggplot(aes(
      x = reorder(k_labelled_dates, order), 
      y = estimate, 
      shape = strata_var,
      colour = comparison, 
      fill = strata_var)) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_line(
      aes(y = line, 
          linetype = strata_var,
          group = line_group), 
      colour = palette,
      alpha = 0.6
    ) +
    geom_linerange(
      aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = position_dodge_val)) +
    geom_point(
      position = position_dodge(width = position_dodge_val)
    ) +
    facet_wrap( ~ subgroup, scales = "free", nrow=2) +
    scale_y_log10(
      name = y_lab_adj,
      breaks = anytest_y1[["breaks"]],
      limits = anytest_y1[["limits"]],
      oob = scales::oob_keep
    ) +
    labs(
      x = x_lab,
      title = "Any SARS-CoV-2 test"
    ) +
    scale_linetype_manual(
      name = NULL,
      values = line_types
    ) +
    scale_shape_manual(
      values = point_shapes
    ) +
    scale_colour_manual(guide = "none",
                        values = palette) +
    scale_fill_manual(guide = "none",
                      values = fill_shapes) +
    guides(
      shape = guide_legend(
        title = NULL, 
        override.aes = list(colour = palette, fill = fill_shapes),
        nrow=leg_rows, byrow=TRUE)
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=10),
      
      axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      strip.text = element_text(size=8),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0, size=10),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = leg_pos,
      legend.key.width = unit(2, 'cm'),
      legend.text = element_text(size=10)
    ) 
  
  # save the plot
  ggsave(plot_vax_anytest,
         filename = here::here(release_folder, glue("hr_vax_anytest_{plot_comparison}_{strata}.png")),
         width=strata_page_width, height=strata_page_height, units="cm")
  
}

for (i in c("BNT162b2", "ChAdOx1", "both")) {
  plot_strata(i, strata = "sex")
  plot_strata(i, strata = "ageband")
}


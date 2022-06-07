################################################################################
# This script:
# - combines estimates from all models into csv for release

################################################################################
library(tidyverse)
library(RColorBrewer)
library(lubridate)
library(glue)

################################################################################
fs::dir_create(here::here("output", "report", "data"))

release_folder <- "release20220226"

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("output", "lib", "study_parameters.rds"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("output", "lib", "outcomes.rds")
)
outcomes_order <- c(3,4,2,5,1)
outcomes_long <- names(outcomes)
outcomes_long[outcomes=="covidadmitted"] <- "COVID-19 hospitalisation"
names(outcomes) <- outcomes_long
rm(outcomes_long)

# read subgroups
subgroups <- readr::read_rds(
  here::here("output", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)
subgroups_order <- c(4,1,3,2)

# define comparisons
comparisons <- c("BNT162b2", "ChAdOx1", "both")

# min and max follow-up dates per subgroup
min_max_fu_dates <- readr::read_rds(
  here::here("output", "lib", glue("data_min_max_fu.rds"))) %>%
  mutate(across(ends_with("date"),
                ~ str_c(day(.x), " ", month(.x, label=TRUE))))

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
K <- study_parameters$max_comparisons
ends <- seq(2, (K+1)*4, 4)
starts <- ends + 1
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")

################################################################################
# read estimates data
estimates_all <- readr::read_csv(
  here::here(release_folder, "estimates_all.csv")) %>%
  mutate(across(model, ~as.integer(str_remove(.x, "unadjusted")))) %>%
  mutate(across(comparison, ~str_replace(.x, "ChAdOx", "ChAdOx1")))

# read metareg data
metareg_results_k <- readr::read_rds(
  here::here(release_folder, "metareg_results_k.rds")) %>%
  select(subgroup, comparison, outcome, k, starts_with("line")) %>%
  mutate(model=2) %>%
  mutate(across(comparison, ~str_replace(.x, "ChAdOx", "ChAdOx1")))

################################################################################

page_height <- 27
page_width <- 16

legend_width <- 15
xlab <- "Weeks since second dose"
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
  # remove as very few events
  filter(!(comparison %in% c("BNT162b2", "both") & subgroup==3 & outcome == "noncoviddeath")) %>%
  mutate(k=as.integer(label)) %>%
  left_join(
    min_max_fu_dates %>%
      mutate(across(subgroup, ~subgroup_labels[subgroups == .x])), 
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
                levels = subgroup_labels[subgroups_order],
                labels = subgroups[subgroups_order]
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
  group_by(k, subgroup, outcome, comparison) %>%
  mutate(order2 = row_number()) %>%
  ungroup() %>%
  mutate(order = order1 + order2) %>%
  left_join(
    metareg_results_k, 
    by = c("subgroup", "comparison", "outcome", "model", "k")
  ) %>%
  mutate(across(model,
                factor,
                levels = 1:2,
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
                levels = subgroups[subgroups_order],
                labels = str_wrap(subgroup_plot_labels[subgroups_order], 25)
  )) %>%
  mutate(line_group = str_c(subgroup, comparison, outcome,  model, sep = "; ")) %>%
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
point_shapes <- 21:24
# breaks and lims for y-axes
primary_vax_y1 <- list(breaks = c(0.02, 0.05, 0.2, 0.5, 1, 2), 
                       limits = c(0.02, 2))
primary_vax_y2 <- list(breaks = c(0,0.5,0.8, 0.95, 0.98))
primary_brand_y1 <- list(breaks = c(0.2, 0.5, 1, 2), 
                         limits = c(0.2, 2))
anytest_y1 <- list(breaks = c(0.5, 1, 2,5,10), 
                   limits = c(0.5, 10))

################################################################################
# vaccine vs unvaccinated
plot_vax <- plot_data %>%
  filter(
    comparison != "both",
    outcome_unlabelled != "anytest",
    as.integer(model) == 2
  ) %>%
  ggplot(aes(
    x = reorder(k_labelled_dates, order), 
    colour = comparison, 
    fill = comparison, 
    shape = subgroup)
  ) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_line(
    aes(y = line, 
        colour=comparison, 
        group = line_group), 
    linetype = "dashed",
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
    name = "Hazard ratio (HR)",
    breaks = primary_vax_y1[["breaks"]],
    limits = primary_vax_y1[["limits"]],
    oob = scales::oob_keep,
    sec.axis = sec_axis(
      ~(1-.),
      name="Estimated vaccine effectiveness = 100*(1-HR)\n ",
      breaks = primary_vax_y2[["breaks"]],
      labels = function(x){formatpercent100(x, 1)}
    )
  ) +
  labs(
    x = str_c(xlab, "**")
  ) +
  scale_colour_discrete(name = NULL) +
  scale_fill_discrete(guide = "none") +
  scale_shape_manual(guide = "none", 
                     values = point_shapes) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=10),
    
    axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 0)),
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
    
    legend.position = c(0.87, 0.14),
    # big margins to cover up grid lines
    legend.margin = margin(t = 30, r = 47, b = 30, l = 40),
    # legend.position = "bottom",
    legend.text = element_text(size=10)
  ) 

# save the plot
ggsave(plot_vax,
       filename = here::here(release_folder, glue("hr_vax.png")),
       width=page_height, height=page_width, units="cm")
ggsave(plot_vax,
       filename = here::here(release_folder, glue("hr_vax.svg")),
       width=page_height, height=page_width, units="cm")

################################################################################
# brand comparison
palette_adj <- gg_color_hue(3, transparency = 1)
i <- 2 # green

plot_brand <- plot_data %>%
  filter(
    comparison == "both",
    outcome_unlabelled != "anytest",
    as.integer(model) == 2
  ) %>%
  droplevels() %>%
  complete(subgroup, comparison, outcome, k_labelled) %>%
  ggplot(aes(
    x = reorder(k_labelled, order), 
    shape = subgroup)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_line(
    aes(y = line, 
        linetype=subgroup, 
        group = line_group), 
    colour = palette_adj[i],
    alpha = 0.6
  ) +
  geom_linerange(
    aes(ymin = conf.low, ymax = conf.high), 
    position = position_dodge(width = position_dodge_val),
    color = palette_adj[i],
    fill = palette_adj[i]
  ) +
  geom_point(
    aes(y = estimate),
    position = position_dodge(width = position_dodge_val),
    color = palette_adj[i],
    fill = palette_adj[i]
  ) +
  facet_grid(outcome ~ ., switch = "y", scales = "free", space = "free_x") +
  scale_y_log10(
    name = "Hazard ratio (HR)",
    breaks = primary_brand_y1[["breaks"]],
    limits = primary_brand_y1[["limits"]],
    oob = scales::oob_keep
  ) +
  labs(
    x = xlab
  ) +
  scale_shape_manual(name = "Subgroup:\n \n HR estimate", values = point_shapes, drop = FALSE) +
  scale_linetype_manual(name="Meta-regression line", values = c("solid", "longdash", "dotted")) +
  guides(shape = guide_legend(order = 1), 
         linetype = guide_legend(order = 2)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=10),
    
    axis.title.x = element_text(size = 10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 0, b = 0, l = 0)),
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
    
    legend.position = "right",
    legend.title = element_text(size=10),
    legend.text = element_text(size=10)
  ) 
ggsave(plot_brand,
       filename = here::here(release_folder, glue("hr_brand.png")),
       width=page_width, height=14, units="cm")
ggsave(plot_brand,
       filename = here::here(release_folder, glue("hr_brand.svg")),
       width=page_width, height=14, units="cm")

################################################################################
# anytest
# vaccine vs unvaccinated
plot_vax_anytest <- plot_data %>%
  filter(
    comparison != "both",
    outcome_unlabelled == "anytest",
    as.integer(model) == 2
  ) %>%
  ggplot(aes(
    x = reorder(k_labelled_dates, order), 
    y = estimate, 
    colour = comparison, 
    fill = comparison, 
    shape = subgroup)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_line(
    aes(y = line, 
        colour=comparison, 
        group = line_group), 
    linetype = "dashed",
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
    name = "Hazard ratio",
    breaks = anytest_y1[["breaks"]],
    limits = anytest_y1[["limits"]],
    oob = scales::oob_keep
  ) +
  labs(
    x = str_c(xlab, "**"),
    title = "Any SARS-CoV-2 test"
  ) +
  scale_colour_discrete(name = NULL) +
  scale_fill_discrete(guide = "none") +
  scale_shape_manual(guide = "none", 
                     values = point_shapes) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=10),
    
    axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 0)),
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
       filename = here::here(release_folder, glue("hr_vax_anytest.png")),
       width=page_width, height=12, units="cm")

################################################################################
# anytest
# brand comparison
palette_adj <- gg_color_hue(3, transparency = 1)
i <- 2 # green

plot_brand_anytest <- plot_data %>%
  filter(
    comparison == "both",
    outcome_unlabelled == "anytest",
    as.integer(model) == 2
  ) %>%
  droplevels() %>%
  complete(subgroup, comparison, outcome, k_labelled) %>%
  ggplot(aes(
    x = reorder(k_labelled, order), 
    y = estimate, 
    shape = subgroup)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_line(
    aes(y = line, 
        linetype=subgroup, 
        group = line_group), 
    colour = palette_adj[i],
    alpha = 0.6
  ) +
  geom_linerange(
    aes(ymin = conf.low, ymax = conf.high), 
    position = position_dodge(width = position_dodge_val),
    color = palette_adj[i],
    fill = palette_adj[i]
  ) +
  geom_point(
    position = position_dodge(width = position_dodge_val),
    color = palette_adj[i],
    fill = palette_adj[i]
  ) +
  scale_y_log10(
    name = "Hazard ratio (HR)",
    breaks = anytest_y1[["breaks"]],
    limits = anytest_y1[["limits"]],
    oob = scales::oob_keep
  ) +
  labs(
    x = xlab,
    title = "Any SARS-CoV-2 test"
  ) +
  scale_shape_manual(name = "Subgroup:\n \n HR estimate", values = point_shapes, drop = FALSE) +
  scale_linetype_manual(name=" \n \n Meta-regression line", values = c("solid", "longdash", "dotted")) +
  guides(shape = guide_legend(order = 1), 
         linetype = guide_legend(order = 2)) +
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
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = "right",
    legend.box = "horizontal",
    legend.title = element_text(size=10),
    legend.text = element_text(size=10)
  ) 
ggsave(plot_brand_anytest,
       filename = here::here(release_folder, glue("hr_brand_anytest.png")),
       width=page_width, height=7, units="cm")

################################################################################
# unadjusted and adjusted estimates for each comparison

plot_unadj_adj <- function(plot_comparison) {
  
  alpha_unadj <- 0.3
  
  # colour palette 
  if (plot_comparison == "both") {
    
    palette_unadj <- gg_color_hue(3, transparency = alpha_unadj)
    palette_adj <- gg_color_hue(3, transparency = 1)
    i <- 2 # green
    
  } else {
    
    palette_unadj <- gg_color_hue(2, transparency = alpha_unadj)
    palette_adj <- gg_color_hue(2, transparency = 1)
    i <- case_when(
      plot_comparison %in% "BNT162b2" ~ 1,  # red 
      plot_comparison %in% "ChAdOx1" ~ 2, # blue
      TRUE ~ NA_real_
    )
    
  }
  
  palette <- c(palette_unadj[i], palette_adj[i])
  
  # vaccine vs unvaccinated
  plot_vax_0 <- plot_data %>%
    filter(
      comparison == plot_comparison,
      outcome_unlabelled != "anytest"
    ) %>%
    ggplot(aes(
      x = reorder(k_labelled_dates, order), 
      y = estimate, 
      colour = model, 
      fill = model, 
      shape = subgroup)) +
    geom_hline(aes(yintercept=1), colour='grey') +
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
        name = "Hazard ratio (HR)",
        breaks = c(0.02, 0.05, 0.2, 0.5, 1, 2),
        limits = c(0.02, 2),
        oob = scales::oob_keep
      )
    
  } else {
    
    plot_vax_1 <- plot_vax_0  +
      scale_y_log10(
        name = "Hazard ratio (HR)",
        breaks = primary_vax_y1[["breaks"]],
        limits = primary_vax_y1[["limits"]],
        oob = scales::oob_keep,
        sec.axis = sec_axis(
          ~(1-.),
          name="Estimated vaccine effectiveness = 100*(1-HR)\n ",
          breaks = c(0,0.5,0.8, 0.95, 0.98),
          labels = function(x){formatpercent100(x, 1)}
        )
      )
    
  }
  
  plot_vax_2 <- plot_vax_1 +
    labs(
      x = str_c(xlab, "**")
    ) +
    scale_colour_manual(name = NULL,
                        values = palette) +
    scale_fill_manual(guide = "none",
                      values = palette) +
    scale_shape_manual(guide = "none", 
                       values = point_shapes) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=10),
      
      axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 0)),
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
         filename = here::here(release_folder, glue("hr_vax_{plot_comparison}.png")),
         width=page_height, height=page_width, units="cm")
  
  # plot comparison
  plot_vax_anytest <- plot_data %>%
    filter(
      comparison == plot_comparison,
      outcome_unlabelled == "anytest"
    ) %>%
    ggplot(aes(
      x = reorder(k_labelled_dates, order), 
      y = estimate, 
      colour = model, 
      fill = model, 
      shape = subgroup)) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_linerange(
      aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = position_dodge_val)) +
    geom_point(
      position = position_dodge(width = position_dodge_val)
    ) +
    facet_wrap( ~ subgroup, scales = "free", nrow=2) +
    scale_y_log10(
      name = "Hazard ratio",
      breaks = anytest_y1[["breaks"]],
      limits = anytest_y1[["limits"]],
      oob = scales::oob_keep
    ) +
    labs(
      x = str_c(xlab, "**"),
      title = "Any SARS-CoV-2 test"
    ) +
    scale_colour_manual(name = NULL,
                        values = palette) +
    scale_fill_manual(guide = "none",
                      values = palette) +
    scale_shape_manual(guide = "none", 
                       values = point_shapes) +
    guides(colour=guide_legend(nrow=2)) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=10),
      
      axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 0)),
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
         filename = here::here(release_folder, glue("hr_vax_anytest_{plot_comparison}.png")),
         width=page_width, height=12, units="cm")
  
}

for (i in c("BNT162b2", "ChAdOx1", "both")) {
  plot_unadj_adj(i)
}


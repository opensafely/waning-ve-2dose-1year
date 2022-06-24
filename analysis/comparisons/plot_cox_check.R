# plot all estimates to check to errors before requesting release of estimates_all.csv

library(tidyverse)
library(glue)

################################################################################
if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")

fs::dir_create(file.path(release_folder, "images"))

################################################################################
cat("-- read subgroups --")
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

################################################################################
cat("-- read estimates_all.csv --")
estimates_all <- readr::read_csv(file.path(release_folder, "estimates_all.csv")) %>%
  filter(variable == "k", label != "0") %>%
  mutate(across(c(estimate, conf.low, conf.high), exp)) %>%
  mutate(across(subgroup, factor, levels = subgroup_labels, labels = subgroups)) 
  
################################################################################
cat("-- define gg_color_hue --")
gg_color_hue <- function(n, transparency = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = transparency)[1:n]
}

################################################################################
cat("-- set parameters for plot --")
position_dodge_val <- 0.8
alpha_unadj <- 0.3
palette_unadj <- c(gg_color_hue(2, transparency = alpha_unadj),
                   gg_color_hue(3, transparency = alpha_unadj)[2])
palette_adj <- c(gg_color_hue(2, transparency = 1),
                 gg_color_hue(3, transparency = 1)[2])

comparisons <- c("BNT162b2", "ChAdOx1", "both")

colour_levs <- c(str_c(comparisons, " unadjusted"), str_c(comparisons, " adjusted"))

palette_all <- c(palette_unadj, palette_adj)
names(palette_all) <- colour_levs

################################################################################
# create function for plot_check
plot_check <- function(c) {
  
  cat(glue("-- comparison = {c} --"))
  
  data_plot <- estimates_all %>%
    filter(
      comparison %in% c,
      outcome != "covidemergency"
      ) %>%
    droplevels()
  
  p <- data_plot %>%
    mutate(
      colourvar = factor(
        glue("{comparison} {model}"),
        levels = colour_levs
        )
      ) %>%
    mutate(across(outcome, 
                  factor,
                  levels = c("covidadmitted", "coviddeath", "postest", "noncoviddeath", "anytest")
                  )) %>%
    mutate(across(label, as.integer)) %>%
    ggplot(aes(
      x = label, 
      y = estimate,
      colour = colourvar
    )) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_linerange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = position_dodge_val)) +
    geom_point(
      position = position_dodge(width = position_dodge_val)
    ) +
    facet_grid(outcome ~ subgroup)  +
    scale_x_continuous(
      breaks = seq(2,12,2)
    ) +
    scale_y_log10(
      name = "HR",
      breaks = c(0.00, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10),
      # limits = c(y_lower, y_upper),
      oob = scales::oob_keep
    ) +
    scale_colour_manual(
      values = palette_all
    ) +
    guides(
      colour = guide_legend(nrow = 2, byrow=TRUE)
    ) +
    labs(
      x = "comaprison period"
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=8),
      
      axis.title.x = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      strip.text = element_text(size=8),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0, size = 8),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      plot.margin = margin(t=10, r=15, b=10, l=10),
      
      legend.position = "bottom"
      
    ) 
  
  cat("-- save plot --")
  # ggsave(p,
  #        filename = here::here("output", "models_cox", "images", glue("plot_check_{c}.svg")),
  #        width=20, height=20, units="cm")
  
  ggsave(p,
         filename = file.path(release_folder, "images", glue("plot_check_{c}.svg")),
         width=20, height=20, units="cm")
}

################################################################################
try(plot_check(c = "BNT162b2"))
try(plot_check(c = "ChAdOx1"))
try(plot_check(c = "both"))
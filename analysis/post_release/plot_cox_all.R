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

if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")

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
min_max_fu_dates <- bind_rows(lapply(
  1:4,
  function(x)
    readr::read_csv(file.path(release_folder, glue("check_fu_{x}.csv"))) %>%
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
estimates_all <- readr::read_csv(here::here(release_folder, "estimates_all.csv")) 

# cumulative incidence data
survtable_redacted <- readr::read_csv(
  here::here(release_folder, "survtable_redacted.csv")) %>%
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
                labels = str_wrap(subgroups, 25)
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
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")

################################################################################

page_height <- 27
page_width <- 16

# y_lab <- "Hazard Ratio (HR)"
# y_lab_2 <- "Estimated vaccine effectiveness = 100*(1-HR)\n "
# y_lab_adj <- "adjusted Hazard Ratio (aHR)"
# y_lab_adj_2 <- "Estimated vaccine effectiveness = 100*(1-aHR)\n "
# x_lab <- "Weeks since second dose"
# # x_lab_nofootnote <- str_remove(x_lab, "\\*\\*")

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
                levels = unname(outcomes[outcomes_order]),
                labels = str_wrap(names(outcomes[outcomes_order]), 10)
  )) %>%
  mutate(across(subgroup,
                factor,
                levels = subgroups,
                labels = str_wrap(subgroup_plot_labels, 25)
  )) %>%
  mutate(line_group = str_c(subgroup, comparison, outcome,  model, sep = "; "))


# spacing of points on plot
position_dodge_val <- 0.6

# colours of points
palette_adj <- RColorBrewer::brewer.pal(3, "Dark2")
names(palette_adj) <- c(comparisons[1:2], "Unvaccinated")

# shapes of points
variants <- c("alpha", "delta", "omicron")
variant_shapes <- c(21,24,22)
names(variant_shapes) <- variants

# fill of points
fill_groups <- c("other_BNT162b2", "other_ChAdOx1", "omicron_BNT162b2", "omicron_ChAdOx1")
palette_fill <- c("white", "white", palette_adj[1:2])
names(palette_fill) <- fill_groups

# breaks and lims for y-axes
primary_vax_y1 <- list(breaks = c(0.02, 0.05, 0.2, 0.5, 1, 2, 4), 
                       limits = c(0.02, 4))
primary_vax_y2 <- list(breaks = c(0,0.5,0.8, 0.95, 0.98))
primary_brand_y1 <- list(breaks = c(0.2, 0.5, 1, 2, 5), 
                         limits = c(0.2, 5))
anytest_y1 <- list(breaks = c(0.5, 1, 2, 5), 
                   limits = c(0.5, 5))

################################################################################
# vaccine vs unvaccinated
plot_vax_data <- plot_data %>%
  filter(
    comparison != "both",
    outcome_unlabelled != "anytest",
    as.integer(model) == 2
  ) %>%
  mutate(across(outcome,
                factor,
                levels = levels(plot_data$outcome),
                labels = str_c("HR for\n", levels(plot_data$outcome)))) %>%
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
  ggplot(aes(
    x = k, 
    colour = comparison, 
    shape = variant,
    # alpha = variant,
    fill = fill_group
  )) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_linerange(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = position_dodge_val)
  ) +
  geom_point(
    aes(y = estimate),
    position = position_dodge(width = position_dodge_val)
  ) +
  facet_grid(
    outcome ~ subgroup, 
    switch = "y",
    # scales = "free",
    space = "free_x"
    ) +
  labs(x = "Weeks since second dose") +
  scale_x_continuous(
    expand = c(0,0),
    breaks = 1:12, # scale is time since start of period 1
    labels = weeks_since_2nd_vax # label scale as time since second vax
  ) +
  scale_y_log10(
    name = NULL,#"HR",
    breaks = primary_vax_y1[["breaks"]],
    limits = primary_vax_y1[["limits"]],
    oob = scales::oob_keep,
    sec.axis = sec_axis(
      ~(1-.),
      name="Estimated VE\n=\n100 x (1-HR)",
      breaks = primary_vax_y2[["breaks"]],
      labels = function(x){formatpercent100(x, 1)}
    )
  ) +
  scale_fill_discrete(guide = "none") +
  scale_shape_manual(
    values = variant_shapes, 
    name = NULL, 
    guide = "none"
  ) +
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
    shape = guide_legend(
      title = "Dominant variant during period:",
      nrow = 1,
      override.aes = list(
        colour = "black",
        fill = c("white", "white", "black")
      )
    ),
    colour = guide_legend(
      title = NULL,#"Vaccination group:",
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        shape = NA,
        size = 3
      )
    )
    ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=8),
    
    axis.title.x = element_text(size=8, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    
    axis.title.y.right = element_text(
      size = 8,
      margin = margin(t = 0, r = 0, b = 0, l = 10),
      angle = 0,
      vjust = 0.5
    ),
    
    axis.text.x = element_text(size=8, angle=90, hjust=1),
    axis.text.y = element_text(size=8),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size=8),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = c(0.88, 0.64),
    legend.spacing.y = unit(0.1, "cm"),
    # big margins to cover up grid lines
    # legend.margin = margin(t = 10, r = 5, b = 10, l = 5),
    legend.key.size = unit(0.3, "cm"),
    # legend.key.width = unit(0.5, 'cm'),
    # legend.position = "bottom",
    legend.title = element_text(size=8),
    legend.text = element_text(size=8),
    legend.box.background = element_rect(colour = "black", fill = "white")
  ) 

# scale for x-axis
x_breaks <- seq(3, 48, 4)
x_labels <- as.character(x_breaks)
alpha_area <- 0.5
# create plot
max_nrisk <- 1000000#max(survtable_redacted$n.risk)

plot_ci_data <- survtable_redacted %>%
  filter(!(as.integer(subgroup) == 4L & arm == "ChAdOx1")) %>%
  filter(time <=50) %>%
  mutate(across(subgroup,
                factor,
                levels = levels(survtable_redacted$subgroup),
                labels = str_replace(str_replace(levels(survtable_redacted$subgroup), "\\n", " "), "and", "&")
                )) %>%
  mutate(
    nrisk_scaled = n.risk/max_nrisk,
    nrisk_scaled_BNT = if_else(arm == "BNT162b2", nrisk_scaled, NA_real_),
    nrisk_scaled_AZ = if_else(arm == "ChAdOx1", nrisk_scaled, NA_real_),
    nrisk_scaled_unvax = if_else(arm == "Unvaccinated", nrisk_scaled, NA_real_),
  ) %>%
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
  ) %>%
  mutate(y = str_wrap("Cumulative incidence of subsequent dose", 14))

plot_ci <- plot_ci_data %>%
  ggplot(aes(
    x = time,
    # y = c.inc,
    colour = arm
    )) +
  geom_area(
    aes(y = nrisk_scaled_unvax),
    fill = palette_adj["Unvaccinated"],
    linetype = "blank",
    alpha = alpha_area
  ) +
  geom_area(
    aes(y = nrisk_scaled_AZ),
    fill = palette_adj["ChAdOx1"],
    linetype = "blank",
    alpha = alpha_area
  ) +
  geom_area(
    aes(y = nrisk_scaled_BNT),
    fill = palette_adj["BNT162b2"],
    linetype = "blank",
    alpha = alpha_area
  ) +
  geom_line(
    aes(y = c.inc_alphadelta),
    linetype = "dashed",
    # size=1
    ) +
  geom_line(
    aes(y = c.inc_omicron),
    # size=1
  ) +
  facet_grid(
    . ~ subgroup,
    switch = "y",
    space = "free_x"
    ) +
  scale_x_continuous(
    expand = c(0,0),
    breaks = seq(2,50,4),
    labels = NULL
  ) +
  scale_y_continuous(
    name = "Millions of\npeople at risk\n(shaded)",
    limits = c(0,1),
    labels = format(seq(0,1,0.25), nsmall=2),
    oob = scales::oob_keep,
    minor_breaks = seq(0,1,0.125),
    sec.axis = sec_axis(
      ~(.),
      name = str_wrap("Cumulative incidence of subsequent dose (line)", 14),
      breaks = seq(0,1,0.25),
      labels = scales::percent_format()
      # labels = function(x){x*max_nrisk/1000000}
    )
  ) +
  scale_colour_manual(
    name = NULL,
    values = palette_adj,
    guide = "none"
  ) +
  labs(x = NULL) +
  theme_bw() +
  theme(
    
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(
      size=8,
      angle = 90
      ),
    
    axis.title.x = element_text(
      size = 10, 
      margin = margin(t = 20, r = 0, b = 10, l = 0)
      ),
    axis.title.y.left = element_text(
      size = 8, 
      margin = margin(t = 0, r = 10, b = 0, l = 0),
      angle = 0,
      vjust = 0.5
      ),
    axis.title.y.right = element_text(
      size = 8,
      margin = margin(t = 0, r = 0, b = 0, l = 10),
      angle = 0,
      vjust = 0.5
    ),
    
    panel.grid.minor.x = element_blank(),
    # panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size=8),
    axis.ticks.x = element_blank(),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic")
    
  )


plot_combined <- plot_grid(
  plot_ci, plot_vax, 
  nrow = 2, rel_heights = c(0.22,0.78),
  labels = c("A", "B"),
  align="v", axis = c("lr")
)

print(caption_str)

# save the plot
ggsave(plot_combined,
       filename = here::here(release_folder, "images", glue("hr_vax_ci.png")),
       width=page_height, height=page_width, units="cm")

# ggsave(plot_vax + theme(plot.margin = margin(2, 2, 2, 2, "cm")),
#        filename = here::here(release_folder, "images", glue("hr_vax.pdf")),
#        width=page_height, height=page_width, units="cm")


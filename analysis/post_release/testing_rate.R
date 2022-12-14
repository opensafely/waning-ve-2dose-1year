library(tidyverse)
library(glue)

# create and define release folder
release_folder <- here::here("release20221006")

## read and derive metadata 
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# labels for comparison periods
K <- study_parameters$K
ends <- seq(2, (K+1)*4, 4)
starts <- ends + 1
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))
outcomes <- outcomes[!(outcomes %in% "covidemergency")]
# outcomes <- outcomes[!(outcomes %in% c("covidemergency", "anytest"))]
old_names <- names(outcomes)
new_names <- str_replace(old_names, "Any", "any")
new_names <- str_replace(new_names, "Positive", "positive")
new_names <- str_replace(new_names, "Non", "non")
new_names <- str_remove(new_names, " \\(APCS\\)")
names(outcomes) <- new_names

outcomes_order <- c(3,4,2,5,1)
# outcomes_order <- c(2,3,1,4)

cat("read event counts data")
event_counts <- readr::read_csv(file.path(release_folder, "event_counts_all.csv")) 

tests <- c("anytest", "postest")

event_counts_processed <- event_counts %>%
  mutate(across(subgroup, as.integer)) %>%
  filter(subgroup==4, str_detect(outcome, "\\w+test")) 
  


estimates_k <- readr::read_csv(file.path(release_folder, "estimates_all.csv")) 

estimates_k_processed <- estimates_k %>%
  filter(variable == "k", !reference_row, model == "adjusted") %>%
  filter(subgroup==4, str_detect(outcome, "\\w+test")) %>%
  select(-variable, -label, -reference_row) %>%
  mutate(across(c("estimate", "conf.low", "conf.high"), exp))


plot_data <- event_counts_processed %>%
  left_join(estimates_k_processed, by = c("subgroup", "arm" = "comparison", "outcome", "k" = "period")) %>%
  mutate(across(outcome, factor, levels = tests, labels = names(outcomes[outcomes%in%tests]))) %>%
  mutate(k_labelled=k) %>%
  mutate(across(k_labelled, factor, levels=seq_along(weeks_since_2nd_vax), labels = weeks_since_2nd_vax)) %>%
  mutate(group = factor(if_else(arm=="unvax", "Unvaccinated", arm))) %>%
  mutate(person28days=365*person_years/28) %>%
  mutate(rate = events/person28days, plot_group=str_c(outcome, group, sep="_")) %>%
  droplevels() %>%
  mutate(
    variant = factor(
      case_when(
        as.integer(subgroup) == 1L & k == 1 ~ "alpha",
        as.integer(subgroup) == 1L & k < 9 ~ "delta",
        as.integer(subgroup) == 1L ~ "omicron",
        as.integer(subgroup) == 2L & k < 8 ~ "delta",
        as.integer(subgroup) == 2L ~ "omicron",
        as.integer(subgroup) == 3L & k < 7 ~ "delta",
        as.integer(subgroup) == 3L ~ "omicron",
        as.integer(subgroup) == 4L & k < 4 ~ "delta",
        as.integer(subgroup) == 4L ~ "omicron",
        TRUE ~ NA_character_
      ),
      levels = variants
    )
  ) 


# shapes of points
variants <- c("alpha", "delta", "omicron")
variant_shapes <- c(16,17,15)
names(variant_shapes) <- variants

# spacing of points on plot
position_dodge_val <- 0.6

# breaks and lims for y-axes
primary_vax_y1 <- list(breaks = c(0.25, 0.5, 1, 2, 4), 
                       limits = c(0.25, 4))


# colours of points
palette_colour <- RColorBrewer::brewer.pal(3, "Dark2")[1]
palette_alpha <- c(0.4, 0.8)
names(palette_alpha) <- levels(plot_data$outcome)

plot_data %>%
  filter(!is.na(estimate)) %>%
  ggplot(aes(x=k, y = estimate, alpha = outcome)) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_linerange(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(width = position_dodge_val), 
    colour = palette_colour
  ) +
  geom_point(
    aes(y = estimate, shape = variant),
    colour = palette_colour,
    size=2,
    position = position_dodge(width = position_dodge_val)
  ) +
  labs(
    x = "Weeks since second dose",
    y = "Adjusted hazard ratio (log scale)",
    title = "18-39 years (and not clinically vulnerable) subgroup"
    ) +
  scale_x_continuous(
    expand = c(0,0),
    breaks = 1:12, # scale is time since start of period 1
    labels = weeks_since_2nd_vax # label scale as time since second vax
  ) +
  scale_y_log10(
    breaks = primary_vax_y1[["breaks"]],
    limits = primary_vax_y1[["limits"]],
    oob = scales::oob_keep
  ) +
  scale_alpha_manual(
    name = "Outcome",
    values = palette_alpha
    ) +
  scale_shape_manual(
    values = variant_shapes, 
    name = NULL, 
    guide = "none"
  ) +
  guides(
    shape = guide_legend(
      title = "Dominant variant during period:",
      nrow = 1#,
      # override.aes = list(
      #   colour = "black",
      #   fill = c("white", "white", "black")
      # )
    ),
    alpha = guide_legend(
      title = "Outcome:",
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        shape = 16#,
        # size = 3
      )
    )
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    
    legend.position = c(0.75, 0.18),
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size=12),
    legend.text = element_text(size=12),
    legend.box.background = element_rect(colour = "black", fill = "white")
    
  )
ggsave(
  filename = "testing_hrs_1839.png",
  path = file.path(release_folder, "images"),
  width = 15, height=15, units="cm"
)


plot_data %>%
  ggplot(aes(x = k,
             y = rate,
             colour = outcome,
             linetype = group,
             group=plot_group)) +
  geom_line() +
  geom_point() +
  labs(x="weeks since second dose", y="rate per 28 person-days") +
  theme_bw() +
  theme(
    axis.title.x = element_text(margin = margin(t=10)),
    axis.title.y = element_text(margin = margin(r=10)),
    legend.box.background = element_rect(colour = "black"),
    legend.position = c(0.8,0.8)
    )
ggsave(
  filename = "testing_rate.png",
  path = file.path(release_folder, "images"),
  width = 15, height=15, units="cm"
)



plot_data %>%
  filter(arm=="BNT162b2", str_detect(outcome, "any")) %>%
  mutate(min=min(estimate), max=max(estimate)) %>%
  filter(estimate==min|estimate==max) %>%
  mutate(across(c(estimate, conf.low, conf.high), round, 2)) %>%
  transmute(value = glue("{format(estimate, nsmall=2)} [{format(conf.low, nsmall=2)}, {format(conf.high, nsmall=2)}]"))


plot_data %>%
  filter(arm=="BNT162b2", str_detect(outcome, "positive"), k%in%c(1,5)) %>%
  mutate(across(c(estimate, conf.low, conf.high), round, 2)) %>%
  transmute(value = glue("{format(estimate, nsmall=2)} [{format(conf.low, nsmall=2)}, {format(conf.high, nsmall=2)}]"))
  




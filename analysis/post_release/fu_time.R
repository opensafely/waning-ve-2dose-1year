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

data_fu <- data_fu %>%
  # check that the subgroup x vaccine groups are correct 
  # if not work out if re-extraction is necessary 
  # (might be possible to work it out from the order??)
  group_by(subgroup, k) %>%
  mutate(group = row_number()) %>%
  ungroup() %>%
  # remove the BNT162b2 group from subgroup 3 
  filter(!(subgroup == 3 & group == 2)) %>%
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

data_fu %>%
  rename(Subgroup=subgroup) %>%
  ggplot(aes(x = date, y = n, colour = Subgroup, fill = Subgroup, group = glue("{Subgroup}, {k}"))) +
  geom_area(alpha=0.25, position = "identity") +
  # geom_vline(xintercept = c(delta_start, omicron_start, delta_end), linetype = unname(linetype_palette)) +
  geom_vline(aes(xintercept = delta_start, linetype = names(linetype_palette)[1])) +
  geom_vline(aes(xintercept = omicron_start, linetype = names(linetype_palette)[2])) +
  geom_vline(aes(xintercept = delta_end, linetype = names(linetype_palette)[3])) +
  facet_grid(rows = "k") +
  scale_x_date(
    expand = c(0,0)
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
    strip.background = element_blank(),
    legend.position = c(0.275,0.2),
    legend.box.background = element_rect(colour = "black", fill = "white")
  )

ggsave(
  filename = file.path(image_path, glue("fu_all.png")),
  width=15, height=20, units="cm"
  )


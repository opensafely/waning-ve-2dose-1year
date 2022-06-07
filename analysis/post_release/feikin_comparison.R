# In our main analysis we have quantified waning in terms of ratios of hazard ratio.
# As previous studies (e.g. Feikin et al.) have done so using the percentage point different in 
# vaccine effectiveness between time-points, we would like to quantify waning 
# in this way to facilitate cross-study comparisons.

library(tidyverse)

estimates_all <- readr::read_csv(here::here("release20220505", "estimates_all.csv"))

# calculate 95% CI for percentage point difference in vaccine effectiveness
data <- estimates_all %>%
  filter(
    variable=="k",
    subgroup == "1",
    period %in% c(1,6),
    comparison != "both",
    outcome %in% c("covidadmitted", "postest"),
    !reference_row,
    model == "adjusted"
  ) %>%
  rename(loghr = estimate) %>%
  mutate(
    seloghr = (conf.high - conf.low)/(2*qnorm(0.975)),
    hr = exp(loghr),
    sehr = hr*seloghr
    ) %>%
  select(comparison, outcome, period, hr, sehr) %>%
  pivot_wider(
    names_from = period,
    values_from = c(hr, sehr)
  ) %>%
  mutate(
    # difference between hazard ratios estimated at 1 and 6 months after 2nd dose
    hr_diff = hr_6 - hr_1,
    # standard error for hr_diff derived using sum of squares
    # is this ok...hr_1 and hr_6 from same sample at difference time points
    sehr_diff = sqrt(sehr_6^2 + sehr_1^2),
    # lower and upper bounds for confidence interval
    hr_diff_lower = hr_diff - sehr_diff*qnorm(0.975),
    hr_diff_upper = hr_diff + sehr_diff*qnorm(0.975),
    # transform to vaccine effectiveness (percent) scale
    ve_diff = 100*hr_diff,
    ve_lower = 100*hr_diff_lower,
    ve_upper = 100*hr_diff_upper
    ) 

# format to paste into paper
data %>%
  mutate(across(c(ve_diff, ve_lower, ve_upper), format, trim=TRUE, digits=2, nsmall=2)) %>%
  transmute(
    comparison,
    outcome,
    value = glue("{ve_diff} ({ve_lower} to {ve_upper})")
  )


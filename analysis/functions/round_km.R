# Function from Will's osutils package

ceiling_any <- function(x, to=1){
  # round to nearest 100 millionth to avoid floating point errors
  ceiling(plyr::round_any(x/to, 1/100000000))*to
}
# Kaplan-Meier estimates, rounded for disclosure control ----
#' Rounded Kaplan-Meier curves
#'
#' @param data A data frame containing the required survival times
#' @param time Event/censoring time variable, supplied as a character. Must be numeric >0
#' @param event Event indicator variables supplied as a character. Censored (`0`/`FALSE`) or not (`1`/`TRUE`). Must be logical or integer with values zero or one
#' @param strata names of stratification / grouping variables, supplied as a character vector of variable names
#' @param threshold Redact threshold to apply
#' @return A tibble with rounded numbers of at risk, events, censored, and derived survival estimates, by strata
#'
#' @details  This function rounds Kaplan-Meier survival function estimates by delaying events times until at least `threshold` events have occurs.
#'
#' @export
#' 
round_km <- function(data, time, event, strata, threshold=6){
  stopifnot("Missing values not allow in `time`" = all(!is.na(data[[time]])))
  stopifnot("Missing values not allow in `event`" = all(!is.na(data[[event]])))
  dat_surv <-
    data %>%
    select(all_of(c(time, event, strata))) %>%
    rename(time={{ time }}, event = {{ event }}) %>%
    group_by(across(all_of(strata))) %>%
    nest() %>%
    mutate(
      n_events = map_int(data, ~sum(.x[[event]], na.rm=TRUE)),
      surv_obj = map(data, ~{
        survfit(Surv(time, event) ~ 1, data = .x, conf.type="log-log")
      }),
      surv_obj_tidy = map(surv_obj, ~broom::tidy(.x)),
    ) %>%
    select(strata, n_events, surv_obj_tidy) %>%
    unnest(surv_obj_tidy)
  dat_surv_rounded <-
    dat_surv %>%
    mutate(
      # Use ceiling not round. This is slightly biased upwards,
      # but means there's no disclosure risk at the boundaries (0 and 1) where masking would otherwise be threshold/2
      surv = ceiling_any(estimate, 1/floor(max(n.risk, na.rm=TRUE)/(threshold))),
      surv.ll = ceiling_any(conf.low, 1/floor(max(n.risk, na.rm=TRUE)/(threshold))),
      surv.ul = ceiling_any(conf.high, 1/floor(max(n.risk, na.rm=TRUE)/(threshold))),
      cml.event = ceiling_any(cumsum(n.event), threshold),
      cml.censor = ceiling_any(cumsum(n.censor), threshold),
      n.event = c(NA, diff(cml.event)),
      n.censor = c(NA, diff(cml.censor)),
      n.risk = ceiling_any(max(n.risk, na.rm=TRUE), threshold) - (cml.event + cml.censor)
    ) %>%
    select(all_of(strata), time, surv, surv.ll, surv.ul, n.risk, n.event, n.censor)
  dat_surv_rounded
}

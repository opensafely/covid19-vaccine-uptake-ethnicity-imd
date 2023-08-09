# add time variables
add_time_vars <- function(.data) {
  .data %>%
    mutate(
      # Create *_time which is the days between elig_date and *_date
      covid_vax_disease_1_time = as.integer(difftime(covid_vax_disease_1_date, elig_date, units = "days")),
      death_time = as.integer(difftime(death_date, elig_date, units = "days")),
      dereg_time = as.integer(difftime(dereg_date, elig_date, units = "days"))
    ) %>%
    # Replace times > 182 days (26 weeks) with NA, as this is the end of follow-up
    mutate(across(ends_with("_time"), ~if_else(.x <= 26*7, .x, NA_integer_)))
}

add_age_jcvi_group <- function(.data) {
  .data %>%
    mutate(
      age_jcvi_group = cut(
        age_jcvi, c(18,seq(25,80,by=5), 120),
        right = FALSE, include.lowest = TRUE
        )
    )
}

# calculate stats for flowcharts
flow_stats_rounded <- function(.data, to) {
  .data %>%
    mutate(
      n = roundmid_any(n, to = to),
      n_exclude = lag(n) - n,
      pct_exclude = n_exclude/lag(n),
      pct_all = n / first(n),
      pct_step = n / lag(n),
    ) %>%
    mutate(across(starts_with("pct_"), ~round(.x, 3)))
}

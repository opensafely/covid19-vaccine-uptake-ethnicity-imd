# the following code only runs when running locally
# don't need to run when running on TPP as dummy data not needed
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  
  library(tidyverse)
  library(rlang)
  library(lubridate)
  
  source(here::here("analysis", "design.R"))
  
  set.seed(34675)
  
  extract <- arrow::read_feather(here::here("output", "extract", "input.feather"))
  
  # translate conditional statements from SQL to R
  # atrisk
  atrisk_group_r <- str_replace_all(atrisk_group, "OR", "|")
  # jcvi
  jcvi_groups_definition <- str_replace_all(jcvi_groups$definition, "AND" , "&")
  jcvi_groups_definition <- str_replace_all(jcvi_groups_definition, "DEFAULT" , "TRUE")
  jcvi_groups_r <- str_c(
    str_c(jcvi_groups_definition, "~\"", jcvi_groups$group, "\""), 
    collapse = "; "
    )
  # elig_date
  elig_dates_definition <- str_replace_all(elig_dates$description, "OR" , "|")
  elig_dates_definition <- str_replace_all(elig_dates_definition, "AND" , "&")
  elig_dates_definition <- str_replace_all(elig_dates_definition, "=" , "==")
  elig_dates_definition <- str_replace_all(elig_dates_definition, ">==" , ">=")
  elig_dates_definition <- str_replace_all(elig_dates_definition, "DEFAULT" , "TRUE")
  elig_dates_r <- str_c(
    str_c(elig_dates_definition, "~\"", elig_dates$date, "\""), 
    collapse = "; "
  )
  
  dummy_data <- extract %>%
    # age_2 same as age_1
    mutate(across(age_2, ~ age_1)) %>%
    # make sure vaccine dates are sensible:
    mutate(
      across(
        covid_vax_disease_2_date,
        ~if_else(
          !is.na(.x),
          covid_vax_disease_1_date + days(round(rnorm(nrow(extract), mean = 84, sd = 7), 0)), 
          .x
        )
      )
    ) %>%
    # derive at risk group
    mutate(across(atrisk_group, ~ !! parse_expr(atrisk_group_r))) %>%
    # derive jcvi group
    mutate(across(jcvi_group, ~ case_when(!!! parse_exprs(jcvi_groups_r)))) %>%
    mutate(across(jcvi_group, ~ factor(.x, levels = levels(extract$jcvi_group)))) %>%
    # derive eligibility date
    mutate(across(elig_date, ~ as.POSIXct(case_when(!!! parse_exprs(elig_dates_r))))) %>%
    # make corrections to dereg_date
    mutate(across(dereg_date, ~ if_else(.x < elig_date, as.POSIXct(NA_character_), .x)))
  
    arrow::write_feather(
      dummy_data,
      here::here("output", "extract", "dummy_data.feather")
      )
    
} else {
  
  # when running on TPP save empty output to keep project.yaml happy
  arrow::write_feather(
    dplyr::tibble(),
    here::here("output", "extract", "dummy_data.feather")
  )
  
}
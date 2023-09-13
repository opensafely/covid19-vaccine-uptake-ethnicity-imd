# Setup ------------------------------------------------------------------------
# load libraries
library(tidyverse)
library(here)
library(rlang)
library(glue)

# load design script
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "processing.R"))

# create output directory
dir_path <- here("output", "extract")

# load study definition extract
extract <- arrow::read_feather(file.path(dir_path, "input.feather"))

# Fix dummy_data  --------------------------------------------------------------
# don't execute if running on TPP
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {
  
  # load custom dummy data
  dummy_data <- arrow::read_feather(file.path(dir_path, "dummy_data.feather"))  
  
  # run checks
  source(here("analysis", "functions", "check_dummy_data.R"))
  check_dummy_data(studydef = extract, custom = dummy_data)
  
  # replace extract with custom dummy data
  extract <- dummy_data
  
  # clean up
  rm(dummy_data)
  
} 

# summarise extracted data
extract %>% my_skim(path = file.path(dir_path, "skim_extract.txt"))

# Process data  ----------------------------------------------------------------
# initial processing
data_processed <- extract %>%
  # because date types are not returned consistently by cohort extractor
  mutate(across(ends_with("_date"), ~ as.Date(.))) %>%
  # define factor levels
  mutate(
    jcvi_group = factor(jcvi_group, levels = jcvi_groups$group),
    sex = factor(sex, levels = sex_levels),
    region = factor(region, levels = regions$region),
    ethnicity = factor(ethnicity, levels = ethnicity_levels),
    imd_Q5 = factor(imd_Q5, levels = imd_Q5_levels),
    # define age on different dates according to JCVI groups
    # otherwise a few patients will end up in a different age group to the
    # rest of their jcvi group
    age_jcvi = if_else(jcvi_group %in% c("10", "11", "12"), age_2, age_1)
  ) %>%
  # reorder so patient_id is first column
  select(patient_id, everything()) 
  

# tidy up
rm(extract)

# Perform checks  --------------------------------------------------------------
# checks atrisk_group, jcvi_groups & elig_date derived correctly in study def
# atrisk_group
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

check_groups <- data_processed %>%
  mutate(
    atrisk_group_r = !! parse_expr(atrisk_group_r),
    jcvi_group_r = case_when(!!! parse_exprs(jcvi_groups_r)),
    elig_date_r = as.Date(case_when(!!! parse_exprs(elig_dates_r)))
  ) %>%
  transmute(
    check_atrisk_group = atrisk_group == atrisk_group_r,
    check_jcvi_group = as.character(jcvi_group) == jcvi_group_r,
    check_elig_date = elig_date == elig_date_r
  ) %>%
  summarise_all(~sum(!.))

if (any(unlist(check_groups) > 0)) {
  print(check_groups)
  stop("Groups derived in study definition do not match those derived in R.")
}

# Apply eligibility criteria  --------------------------------------------------
data_processed <- data_processed %>%
  mutate(
    # create variables for applying the eligibility criteria
    aged_over_18_grp_12 = (age_2 >= 18) & (jcvi_group %in% "12"),
    aged_over_18_grps_04b_06 = (age_1 >= 18) & (jcvi_group %in% c("04b", "6")),
    alive_on_elig_date = is.na(death_date) | elig_date <= death_date,
    aged_under_120 = age_2 < 120,
    no_vax_before_start = is.na(covid_vax_disease_1_date) | as.Date(study_parameters$start_date) <= covid_vax_disease_1_date,
    sex_recorded = !is.na(sex),
    region_recorded = !is.na(region),
    imd_recorded = !(is.na(imd_Q5) | (imd_Q5 %in% "Unknown")),
    ethnicity_recorded = !(is.na(ethnicity) | (ethnicity %in% "Unknown")),
    # apply the eligibility criteria
    c0_descr = "All patients in OpenSAFELY-TPP",
    c0 = TRUE,
    c1_descr = glue("   aged 18 years or over by {study_parameters$ref_age_2}"),
    c1 = c0 & (aged_over_18_grp_12 | aged_over_18_grps_04b_06),
    c2_descr = "   alive at start of eligibility date",
    c2 = c1 & alive_on_elig_date,
    c3_descr = "   registered with one TPP general practice between 2020-01-01 and start of eligibility date",
    c3 = c2 & has_follow_up,
    c4_descr = glue("   aged under 120 years on {study_parameters$ref_age_2}"),
    c4 = c3 & aged_under_120,
    c5_descr = "   no vaccination before rollout",
    c5 = c4 & no_vax_before_start,
    c6_descr = "   sex, region, IMD and ethnicity recorded at start of eligibility date",
    c6 = c5 & sex_recorded & region_recorded & imd_recorded & ethnicity_recorded,
    include = c6
  )

# data for flowchart
data_flowchart <- data_processed %>%
  select(patient_id, matches("^c\\d+")) %>%
  rename_at(vars(matches("^c\\d+$")), ~str_c(., "_value")) %>%
  pivot_longer(
    cols = matches("^c\\d+"),
    names_to = c("crit", ".value"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  group_by(crit, descr) %>%
  summarise(n = sum(value), .groups = "keep") %>%
  ungroup() %>%
  rename(criteria = descr) %>%
  arrange(crit) 

# save flowchart without rounding in case needed for debugging
data_flowchart %>%
  flow_stats_rounded(1) %>%
  write_csv(file.path(dir_path, "flowchart_raw.csv"))

# save flowchart with rounding for releasing
data_flowchart %>%
  flow_stats_rounded(threshold) %>%
  write_csv(file.path(dir_path, glue("flowchart_midpoint{threshold}.csv")))

# only keep the variables that are needed for the eligible patients 
data_eligible <- data_processed %>%
  filter(include) %>%
  select(
    patient_id, elig_date, covid_vax_disease_1_date, death_date, dereg_date, 
    age_jcvi, jcvi_group, region, sex, imd_Q5, ethnicity
    ) 

# tidy up
rm(data_processed)

# summarise eligible data
data_eligible %>% my_skim(path = file.path(dir_path, "skim_eligible.txt"))

# save data_eligible
write_rds(
  data_eligible,
  file.path(dir_path, "data_eligible.rds"),
  compress = "gz"
  )

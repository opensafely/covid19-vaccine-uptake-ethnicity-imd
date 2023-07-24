# load libraries
library(tidyverse)
library(here)

# load design script
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))

# create output directory
dir_path <- here("output", "extract")

# load study definition extract
extract <- arrow::read_feather(file.path(dir_path, "input.feather"))

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
    imd_Q5 = factor(imd_Q5, levels = imd_Q5_levels)
  ) %>%
  # reorder so patient_id is first column
  select(patient_id, everything()) 
  

# tidy up
rm(extract)

# TODO Elsie
# write in checks to make sure at risk group, JCVI groups and elig dates have 
# been derived correctly in study definition

# apply exclusions
data_processed <- data_processed %>%
  mutate(
    # create variables for applying the eligibility criteria
    aged_over_18 = age_2 >=18,
    alive_on_elig_date = is.na(death_date) | elig_date <= death_date,
    aged_under_120 = age_2 < 120,
    no_vax_before_start = is.na(covid_vax_disease_1_date) | as.Date(study_parameters$start_date) <= covid_vax_disease_1_date,
    sex_recorded = !is.na(sex),
    region_recorded = !is.na(region),
    imd_recorded = imd_Q5 != "Unknown",
    ethnicity_recorded = ethnicity != "Unknown",
    # apply the eligibility criteria
    c0_descr = "All patients in OpenSAFELY-TPP",
    c0 = TRUE,
    c1_descr = "   aged 18 years or over",
    c1 = c0 & aged_over_18,
    c2_descr = "   alive on eligibility date",
    c2 = c1 & alive_on_elig_date,
    c3_descr = "   registered with one TPP general practice between 2020-01-01 and eligibility date",
    c3 = c2 & has_follow_up,
    c4_descr = "   aged under 120 years",
    c4 = c3 & aged_under_120,
    c5_descr = "   no vaccination before rollout",
    c5 = c4 & no_vax_before_start,
    c6_descr = "   sex, region, IMD and ethnicity recorded",
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
    age_1, jcvi_group, region, sex, imd_Q5, ethnicity
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

################################################################################
# This script:
# creates metadata for aspects of the study design
################################################################################

# Import libraries ----
library(tidyverse)
library(lubridate)
library(glue)

# create lib directory ----
fs::dir_create(here::here("analysis", "lib"))

# threshold for midpoint rounding
threshold <- 6

################################################################################
# create study_parameters ----
study_parameters <-
  lst(
    seed = 123456L,
    ref_age_1 = "2021-03-31", # reference date for calculating age for phase 1 groups
    ref_age_2 = "2021-07-01", # reference date for calculating age for phase 2 groups
    ref_cev = "2021-01-18", # reference date for calculating eligibility for phase 1 group 4 (CEV)
    ref_ar = "2021-02-15", # reference date for calculating eligibility for phase 1 group 5 (at-risk)
    pandemic_start = "2020-01-01", # rough start date for pandemic in UK
    start_date = "2020-12-08", # start of phase 1 vaccinations
  ) 

study_parameters %>% 
  jsonlite::write_json(
    path = here::here("analysis", "lib", "study_parameters.json"),
    auto_unbox = TRUE, pretty=TRUE
    )

################################################################################
# at risk group definition ----
atrisk_group <- "immuno_group OR
                 ckd_group OR
                 resp_group OR
                 asthma_group OR
                 diab_group OR
                 cld_group OR
                 cns_group OR
                 chd_group OR
                 spln_group OR
                 learndis_group OR
                 sevment_group OR
                 sevobese_group" 
atrisk_group <- str_replace_all(str_remove_all(atrisk_group, "\\n"), "\\s+", " ")

tibble(atrisk_group = atrisk_group) %>%
  readr::write_csv(here::here("analysis", "lib", "atrisk_group.csv"))

################################################################################
# create jcvi_groups ----
jcvi_groups <- 
tribble(
    ~group, ~definition,
    "01", "longres_group AND age_1 > 65",
    "02", "age_1 >=80",
    "03", "age_1 >=75",
    "04a", "age_1 >=70",
    "04b", "cev_group AND age_1 >=16", # TODO change to age_1 >= 18
    "05", "age_1 >=65",
    "06", "atrisk_group AND age_1 >=16", # TODO change to age_1 >= 18
    "07", "age_1 >=60",
    "08", "age_1 >=55",
    "09", "age_1 >=50",
    "10", "age_2 >=40",
    "11", "age_2 >=30",
    "12", "age_2 >=18",
    "99", "DEFAULT",
)

readr::write_csv(jcvi_groups, here::here("analysis", "lib", "jcvi_groups.csv"))

################################################################################
# create elig_dates ----
# group elig_date if within 7 days of previous elig_date (within jcvi_group)
elig_dates <-
tribble(
    ~date, ~description, ~jcvi_groups,
    "2020-12-08", "jcvi_group='01' OR jcvi_group='02'", "01, 02", 
    "2021-01-18", "jcvi_group='03' OR jcvi_group='04a' OR jcvi_group='04b'", "03, 04a, 04b",
    ###
    "2021-02-15", "jcvi_group='05' OR jcvi_group='06'", "05, 06",
    ###
    "2021-02-22", "age_1 >= 64 AND age_1 < 65", "07", 
    "2021-03-01", "age_1 >= 60 AND age_1 < 64", "07",
    ###
    "2021-03-08", "age_1 >= 56 AND age_1 < 60", "08",
    "2021-03-09", "age_1 >= 55 AND age_1 < 56", "08",
    ###
    "2021-03-19", "age_1 >= 50 AND age_1 < 55", "09",
    ###
    "2021-04-13", "age_2 >= 45 AND age_1 < 50", "10",
    "2021-04-26", "age_2 >= 44 AND age_1 < 45", "10",
    "2021-04-27", "age_2 >= 42 AND age_1 < 44", "10",
    "2021-04-30", "age_2 >= 40 AND age_1 < 42", "10",
    ###
    "2021-05-13", "age_2 >= 38 AND age_2 < 40", "11",
    "2021-05-19", "age_2 >= 36 AND age_2 < 38", "11",
    "2021-05-21", "age_2 >= 34 AND age_2 < 36", "11",
    "2021-05-25", "age_2 >= 32 AND age_2 < 34", "11",
    "2021-05-26", "age_2 >= 30 AND age_2 < 32", "11",
    ###
    "2021-06-08", "age_2 >= 25 AND age_2 < 30", "12",
    "2021-06-15", "age_2 >= 23 AND age_2 < 25", "12",
    "2021-06-16", "age_2 >= 21 AND age_2 < 23", "12",
    "2021-06-18", "age_2 >= 18 AND age_2 < 21", "12",
    "2100-12-31", "DEFAULT", "DEFAULT",
) 

readr::write_csv(elig_dates, here::here("analysis", "lib", "elig_dates.csv"))

################################################################################
# create regions ----
regions <- tribble(
  ~region, ~ratio,
  "North East", 0.1,
  "North West", 0.1,
  "Yorkshire and The Humber", 0.1,
  "East Midlands", 0.1,
  "West Midlands", 0.1,
  "East", 0.1,
  "London", 0.2,
  "South West", 0.1,
  "South East", 0.1
)

readr::write_csv(regions, here::here("analysis", "lib", "regions.csv"))

################################################################################
sex_levels <- c("F", "M")
ethnicity_levels <- c("White", "Black or Black British", "Asian or Asian British", "Mixed", "Other", "Unknown")
imd_Q5_levels <- c("1 (most deprived)", "2", "3", "4", "5 (least deprived)")

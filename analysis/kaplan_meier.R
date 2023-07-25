# the following code does not calculate Kaplan-Meier estimates, 
# but I've put it here for now as it will be useful for the groupings

# Count how many patients there are in the imd subgroups by ethnicity ----------
data_imd_eth <- data_eligible %>%
  group_by(imd_Q5, ethnicity) %>%
  summarise(n = n(), .groups = "drop")

data_imd_eth_sex <- data_eligible %>%
  group_by(imd_Q5, ethnicity, sex) %>%
  summarise(n = n(), .groups = "drop")

data_imd_eth_reg <- data_eligible %>%
  group_by(imd_Q5, ethnicity, region) %>%
  summarise(n = n(), .groups = "drop")

data_imd_eth_jcvi <- data_eligible %>%
  group_by(imd_Q5, ethnicity, jcvi_group) %>%
  summarise(n = n(), .groups = "drop")

# Bind these together and save as csv file so we can release from opensafely
data_counts <- bind_rows(
  data_imd_eth %>% mutate(variable = NA_character_, level = NA_character_),
  data_imd_eth_sex %>% mutate(variable = "sex") %>% rename(level = sex),
  data_imd_eth_reg %>% mutate(variable = "region") %>% rename(level = region),
  data_imd_eth_jcvi %>% mutate(variable = "jcvi_group") %>% rename(level = jcvi_group)
) %>%
  # we round the counts using midpoint rounding to reduce risk of secondary disclosure
  # threshold is defined in design.R
  mutate(across(n, ~roundmid_any(.x, to = threshold)))

# Save to .csv file for release
readr::write_csv(
  data_counts, 
  file.path(outdir, glue("group_counts_midpoint{threshold}.csv"))
)

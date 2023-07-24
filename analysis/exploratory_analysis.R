# Setup ------------------------------------------------------------------------

# Load libraries 
library(lubridate)
library(purrr)
library(gtsummary)
library(ggplot2)
library(viridis)
library(here)
library(fs)

# Create output directory
outdir <- here("output", "exploratory")
fs::dir_create(outdir)

# Source metadata and functions
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))

# Import the data 
data_eligible <- readRDS(here("output", "extract", "data_eligible.rds"))

# Summary of the data
# summary(data_eligible)

# Process data ----------------------------------------------------------------
data_eligible <- data_eligible %>%
  mutate(
    # Create *_time which is the days between elig_date and *_date
    covid_vax_disease_1_time = as.integer(difftime(covid_vax_disease_1_date, elig_date, units = "days")),
    death_time = as.integer(difftime(death_date, elig_date, units = "days")),
    dereg_time = as.integer(difftime(dereg_date, elig_date, units = "days"))
    ) %>%
  # replace times > 182 days (26 weeks) with NA, as this is the end of follow-up
  mutate(across(ends_with("_time"), ~if_else(.x <= 26*7, .x, NA_integer_)))

# Check for negative times
cat("\nCheck for negative values (should be 0):\n")
cat("\n`death_time`:\n")
sum(data_eligible$death_time < 0, na.rm = TRUE) # should be 0
cat("\n`dereg_time`:\n")
sum(data_eligible$dereg_time < 0, na.rm = TRUE) # should be 0

# Derive a variable that corresponds to the rank of elig_date within jcvi_group
data_elig_date_rank <- data_eligible %>%
  distinct(jcvi_group, elig_date) %>%
  arrange(jcvi_group, elig_date) %>%
  group_by(jcvi_group) %>%
  mutate(elig_date_rank = factor(rank(elig_date))) %>%
  ungroup()

# do the join here so we only have to do it once, becasue joins are slow
data_eligible <- data_eligible %>%
  left_join(data_elig_date_rank, by = c("elig_date", "jcvi_group"))

rm(data_elig_date_rank)

# Explore distribution of covid_vax_disease_1_time -----------------------------

# Plot using geom_freqpoly()
data_eligible %>%
  ggplot(aes(x = covid_vax_disease_1_time, y = after_stat(count), color = elig_date_rank)) +
  geom_freqpoly(binwidth = 1) +
  facet_wrap(~jcvi_group, scales = "free_y", ncol=4) +
  scale_color_viridis_d() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Distribution of covid_vax_disease_1_time across different eligibility dates",
       x = "Days between eligibility and first vaccination",
       y = "Frequency",
       color = "Eligibility Date")
# Save
ggsave(file.path(outdir, "vax_dates_freqpoly.png"))

# Create the data for replicating this plot locally
# Counts of individuals vaccinated on each day, grouped by jcvi_group and elig_date
data_vax_counts <- data_eligible %>%
  # get rid of individuals who did not get vaccinated during follow-up
  filter(!is.na(covid_vax_disease_1_time)) %>%
  group_by(jcvi_group, elig_date, elig_date_rank, covid_vax_disease_1_time) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(n = roundmid_any(n, to = threshold))  # Apply rounding to the counts

# head(data_vax_counts)

# Save to .csv file for release
readr::write_csv(
  data_vax_counts,
  file.path(outdir, glue("data_vax_counts_midpoint{threshold}.csv"))
  )

# Plot using geom_line()
data_vax_counts %>%
  ggplot(aes(x = covid_vax_disease_1_time, y = n, color = elig_date_rank)) +
  geom_line() +
  facet_wrap(~jcvi_group, scales = "free_y", nrow = 4) +
  scale_color_viridis_d() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Distribution of covid_vax_disease_1_date across different eligibility dates",
       x = "Days between eligibility and first vaccination",
       y = "Frequency",
       color = "Eligibility Date")
# Save
ggsave(file.path(outdir, "vax_dates_line.png"))


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

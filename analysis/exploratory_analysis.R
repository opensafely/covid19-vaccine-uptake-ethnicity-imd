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
source(here("analysis", "functions", "processing.R"))

# Import the data 
data_eligible <- readRDS(here("output", "extract", "data_eligible.rds"))

# Summary of the data
# summary(data_eligible)

# Process data ----------------------------------------------------------------
data_eligible <- data_eligible %>%
  # add *_time vars - see analysis/functions/processing.R for add_time_vars()
  add_time_vars() %>%
  # add age_jcvi_group
  add_age_jcvi_group()

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
# (stratified by jcvi_group and elig_date)

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


# Explore distribution of variables stratified by ethnicity x imd --------------
# variables = age_jcvi_group, sex, region and jcvi_group

# TODO write a function based on the code below that takes a variable as input
# (e.g. age_jcvi_group) and generates data_bar_plot
data_bar_plot <- data_eligible %>%
  group_by(ethnicity, imd_Q5, age_jcvi_group) %>%
  summarise(n = roundmid_any(n(), to = threshold)) %>%
  ungroup(age_jcvi_group) %>%
  mutate(percent = 100*n/sum(n)) %>%
  ungroup()

# TODO write a function based on the code below to create a bar plot 
# (feel free to play around with the appearance)
data_bar_plot %>%
  ggplot(aes(x = age_jcvi_group, y = percent)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(ethnicity), cols = vars(imd_Q5)) 

# TODO save a single dataset that contains data_bar_plot for each variable 
# (you'll have to add a column to indicate which variable the data refers to)



# TODO
# Function for visualization data

generate_data_bar_plot <- function(variable) {
  data_bar_plot <- data_eligible %>%
    group_by(ethnicity, imd_Q5, !!sym(variable)) %>%
    summarise(n = roundmid_any(n(), to = threshold)) %>%
    ungroup(!!sym(variable)) %>%
    mutate(percent = 100*n/sum(n)) %>%
    ungroup()
  
  return(data_bar_plot)
}


# Function for bar plot
#create_bar_plot <- function(data, variable) {
#  ggplot(data, aes(x = !!sym(variable), y = percent)) +
#    geom_bar(stat = "identity") +
#    facet_grid(rows = vars(ethnicity), cols = vars(imd_Q5)) +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1))
#}

# Function for bar plot
create_bar_plot <- function(data, variable) {
  ggplot(data, aes_string(x = variable, y = "percent", fill = variable)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(rows = vars(ethnicity), cols = vars(imd_Q5)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) +
    labs(title = paste("Distribution of", variable, "stratified by ethnicity and IMD"),
         x = variable,
         y = "Percentage (%)",
         fill = variable)
}



# Generating the data 

variables <- c("age_jcvi_group", "sex", "region", "jcvi_group")
data_bar_plots <- lapply(variables, generate_data_bar_plot)

# Generate plots for each data_bar_plot
plots <- lapply(seq_along(variables), function(i) {
  create_bar_plot(data_bar_plots[[i]], variables[i])
})

# View outputs of the plots 
#print(plots[[4]])


# Combine all data_bar_plot into a single table 
combined_data_bar_plot <- dplyr::bind_rows(data_bar_plots)

















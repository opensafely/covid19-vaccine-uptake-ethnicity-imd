# code for the main analysis using Kaplan-Meier estimates ----------------------

# Load libraries 
library(ggplot2)
library(here)
library(survival)
library(survminer)

# Create output directory
outdir <- here("output", "exploratory")
fs::dir_create(outdir)

# Source metadata and functions
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "processing.R"))

# Import the data 
data_eligible <- readRDS(here("output", "extract", "data_eligible.rds"))

# Create a numeric variable to measure time from eligibility date to vaccination
data_eligible$days_to_vaccine <- as.numeric(data_eligible$covid_vax_disease_1_date - data_eligible$elig_date)

# Create a categorical identifier of event (vaccine uptake)
data_eligible$vaccine_uptake <- ifelse(is.na(data_eligible$covid_vax_disease_1_date), 0, 1)

# Create a Survival object for kaplan-meier analysis
surv_obj <- Surv(time = data_eligible$days_to_vaccine, event = data_eligible$vaccine_uptake)

# Fit a survival curve for each ethnicity subgroup
fit_ethnicity <- survfit(surv_obj ~ ethnicity, data = data_eligible)

# Plot the survival curves for each ethnicity subgroup
ggsurvplot(fit_ethnicity, data = data_eligible, risk.table = TRUE)


### IMD sub-groups 
##########

# Fit survival curve for each IMD subgroup
fit_imd <- survfit(surv_obj ~ imd_Q5, data = data_eligible)

# Plot survival curves for each IMD subgroup
ggsurvplot(fit_imd, data = data_eligible, risk.table = TRUE)


### Add censoring to the survival plots
############

# Create a 'status' variable where 1 indicates the event (vaccination) occurred, and 0 indicates censoring
data_eligible$status <- ifelse(is.na(data_eligible$covid_vax_disease_1_date) | 
                                 data_eligible$death_date < data_eligible$covid_vax_disease_1_date | 
                                 data_eligible$dereg_date < data_eligible$covid_vax_disease_1_date, 0, 1)

# Create a 'time' variable that represents the time to event or censoring
data_eligible$time <- with(data_eligible, pmin(covid_vax_disease_1_date, death_date, dereg_date, na.rm = TRUE) - elig_date)
data_eligible$time <- as.numeric(data_eligible$time)

# Create a Surv object with the time and status variables
surv_obj <- Surv(time = data_eligible$time, event = data_eligible$status)

# Fit the Kaplan-Meier survival model
fit <- survfit(surv_obj ~ 1)

# Plot the survival curve
ggsurvplot(fit, data = data_eligible, risk.table = TRUE)


#############
#############

# List of unique ethnicities
ethnicities <- unique(data_eligible$ethnicity)

# Create a plot for each ethnicity
for (ethnicity in ethnicities) {
  # Subset the data for the current ethnicity
  data_subset <- data_eligible[data_eligible$ethnicity == ethnicity, ]
  
  # Create a Surv object for kaplan-meier analysis
  surv_obj <- Surv(time = data_subset$days_to_vaccine, event = data_subset$vaccine_uptake)
  
  # Fit a survival curve
  fit <- survfit(surv_obj ~ 1, data = data_subset)
  
  # Plot the survival curve and save it to a variable
  plot <- ggsurvplot(fit, data = data_subset, risk.table = TRUE,
                     title = paste("Kaplan-Meier plot for", ethnicity))
  
  # Save the plot to a file
  ggsave(paste0("km_plot_", ethnicity, ".png"), plot$plot)
}


































# the following code does not calculate Kaplan-Meier estimates, 
# but I've put it here for now as it will be useful for the groupings

# Count how many patients there are in the imd subgroups by ethnicity ----------
#data_imd_eth <- data_eligible %>%
#  group_by(imd_Q5, ethnicity) %>%
#  summarise(n = n(), .groups = "drop")

#data_imd_eth_sex <- data_eligible %>%
#  group_by(imd_Q5, ethnicity, sex) %>%
#  summarise(n = n(), .groups = "drop")

#data_imd_eth_reg <- data_eligible %>%
#  group_by(imd_Q5, ethnicity, region) %>%
#  summarise(n = n(), .groups = "drop")

#data_imd_eth_jcvi <- data_eligible %>%
#  group_by(imd_Q5, ethnicity, jcvi_group) %>%
#  summarise(n = n(), .groups = "drop")

# Bind these together and save as csv file so we can release from opensafely
#data_counts <- bind_rows(
#  data_imd_eth %>% mutate(variable = NA_character_, level = NA_character_),
#  data_imd_eth_sex %>% mutate(variable = "sex") %>% rename(level = sex),
#  data_imd_eth_reg %>% mutate(variable = "region") %>% rename(level = region),
#  data_imd_eth_jcvi %>% mutate(variable = "jcvi_group") %>% rename(level = jcvi_group)
#) %>%
  # we round the counts using midpoint rounding to reduce risk of secondary disclosure
  # threshold is defined in design.R
#  mutate(across(n, ~roundmid_any(.x, to = threshold)))

# Save to .csv file for release
#readr::write_csv(
 # data_counts, 
  #file.path(outdir, glue("group_counts_midpoint{threshold}.csv"))
#)

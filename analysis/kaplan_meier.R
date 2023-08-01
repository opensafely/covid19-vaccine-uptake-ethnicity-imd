# code for the main analysis using Kaplan-Meier estimates ----------------------

# Load libraries 
library(ggplot2)
library(here)
library(survival)
library(survminer)

# Create output directory
outdir <- here("output", "kaplan_meier")
fs::dir_create(outdir)

# Source metadata and functions
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "processing.R"))

# Import the data 
data_eligible <- readRDS(here("output", "extract", "data_eligible.rds"))

# Process data ----------------------------------------------------------------
data_eligible <- data_eligible %>%
  # add *_time vars - see analysis/functions/processing.R for add_time_vars()
  add_time_vars() %>%
  # add age_jcvi_group
  add_age_jcvi_group() %>%
  # define the variables for using in the survfit formula:
  mutate(
    # combine the ethnicity and IMD categories - we will calculate kaplan-meier
    # estimates in each of these categories - you can do this by adding this 
    # variable to the right-hand side of the formula in survfit
    ethnicity_imd = paste0(ethnicity, ", ", imd_Q5),
    # create time to event (tte) variable
    tte = pmin(
      # event time - replace negative times with zero
      pmax(0, covid_vax_disease_1_time), 
      # censoring times
      death_time, dereg_time, 
      # administrative censoring at end of study
      26*7, 
      na.rm = TRUE
      ),
    # status is 1 if the tte corresponds to a vaccination time
    status = covid_vax_disease_1_time == tte
  )

# TODO update the code from here, using the variables that I've defined above.

# Create a Survival object for kaplan-meier analysis
surv_obj <- Surv(time = data_eligible$tte, event = data_eligible$status)

# Fit the Kaplan-Meier survival model
fit <- survfit(surv_obj ~ 1)

# Plot the survival curve
ggsurvplot(fit, data = data_eligible, risk.table = TRUE)

#########
#########


# Fit the Kaplan-Meier survival model for each ethnicity and IMD subgroup
fit_ethnicity_imd <- survfit(surv_obj ~ ethnicity_imd, data = data_eligible)

# Plot the survival curves for each ethnicity and IMD subgroup
ggsurvplot(fit_ethnicity_imd, data = data_eligible, risk.table = TRUE)


##########
##########

# List of unique ethnicity and IMD subgroups
ethnicity_imd_groups <- unique(data_eligible$ethnicity_imd)

# Create a plot for each ethnicity and IMD subgroup
for (group in ethnicity_imd_groups) {
  # Subset the data for the current group
  data_subset <- data_eligible[data_eligible$ethnicity_imd == group, ]
  
  # Create a Surv object for kaplan-meier analysis
  surv_obj <- Surv(time = data_subset$tte, event = data_subset$status)
  
  # Fit a survival curve
  fit <- survfit(surv_obj ~ 1, data = data_subset)
  
  # Plot the survival curve and save it to a variable
  plot <- ggsurvplot(fit, data = data_subset, risk.table = TRUE,
                     title = paste("Kaplan-Meier plot for", group),
                     xlab = "Time (days)", ylab = "Survival probability") +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  # Save the plot to a file
  ggsave(paste0("km_plot_", gsub(" ", "_", gsub(",", "", group)), ".png"), plot$plot)
}

########
########

# List of unique ethnicity and IMD subgroups
ethnicity_imd_groups <- unique(data_eligible$ethnicity_imd)

# Create a plot for each ethnicity and IMD subgroup
for (group in ethnicity_imd_groups) {
  # Subset the data for the current group
  data_subset <- data_eligible[data_eligible$ethnicity_imd == group, ]
  
  # Create a Surv object for kaplan-meier analysis
  surv_obj <- Surv(time = data_subset$tte, event = data_subset$status)
  
  # Fit a survival curve
  fit <- survfit(surv_obj ~ 1, data = data_subset)
  
  # Determine the number of levels in the status variable
  num_status_levels <- length(unique(data_subset$status))
  
  # Define the legend labels based on the number of status levels
  if (num_status_levels == 1) {
    legend_labs <- c("Event")
  } else if (num_status_levels == 2) {
    legend_labs <- c("Censored", "Event")
  } else {
    stop("Unexpected number of status levels")
  }
  
  # Plot the survival curve and save it to a variable
  plot <- ggsurvplot(fit, data = data_subset, risk.table = TRUE,
                     title = paste("Kaplan-Meier plot for", group),
                     xlab = "Time (days)", ylab = "Survival probability",
                     legend.title = "Status", legend.labs = legend_labs,
                     palette = c("#E7B800", "#2E9FDF"),
                     theme = theme_classic2())
  
  # Save the plot to a file
  ggsave(paste0("km_plot_", gsub(" ", "_", gsub(",", "", group)), ".png"), plot$plot)
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

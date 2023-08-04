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
data_surv <- data_eligible %>%
  # add *_time vars - see analysis/functions/processing.R for add_time_vars()
  add_time_vars() %>%
  # add age_jcvi_group
  add_age_jcvi_group() %>%
  # define the variables for using in the survfit formula:
  transmute(
    patient_id,
    # keep all the variables that we'll use to define the subgroups
    ethnicity, imd_Q5, region, jcvi_group, sex,
    # time until an outcome event (replace negative times with zero)
    # do not use na.rm=TRUE here, as we want to keep NAs for unvaccinated people
    tte_outcome = pmax(0, covid_vax_disease_1_time),
    # time until a censoring event (administrative cenosring at 26*7)
    # do use na.rm=TRUE here as we want people with missing death and dereg to 
    # to have tte_censor=26*7
    tte_censor = pmin(death_time, dereg_time, 26*7, na.rm = TRUE),
    # tte is the earlier of tte_outcome and tte_censor
    tte = pmin(tte_outcome, tte_censor, na.rm = TRUE),
    # status = 1 if the person was vaccinated before they were censored
    status = as.integer(!is.na(tte_outcome) & (tte_outcome <= tte_censor))
  )

# Create a Survival object for kaplan-meier analysis
surv_obj <- Surv(time = data_surv$tte, event = data_surv$status)

# Fit the Kaplan-Meier survival model for each ethnicity and IMD subgroup
fit_ethnicity_imd <- survfit(surv_obj ~ ethnicity + imd_Q5, data = data_surv)

# Plot the survival curves for each ethnicity and IMD subgroup
ggsurv_obj <- ggsurvplot(fit_ethnicity_imd, data = data_surv, risk.table = TRUE)

# Extract the plot from ggsurv_obj and get rid of the legend using `+ ggplot2::theme()`, then save the plot
ggsurv_obj$plot <- ggsurv_obj$plot + theme(legend.position = "none")
ggsave(file.path(outdir, "km_plot_all_groups.png"), ggsurv_obj$plot)

# Create a new dataset for vaccine coverage

# Create the survival table
surv_table <- ggsurv_obj$data.survplot
# Filter the survival table to only keep time for 12 weeks and 26 weeks
data_coverage <- surv_table %>%
  filter(time %in% c(12*7, 26*7)) %>%
  mutate(coverage = 1 - surv) %>%
  select(strata, time, coverage, n.event, n.censor, std.err)
# Write the data to a CSV file
write_csv(
  data_coverage,
  file.path(outdir, glue::glue("vaccine_coverage_midpoint{threshold}.csv"))
)


##########
# Examining other covariates
##########

# Create a list of covariates
additional_covariates <- c("region", "jcvi_group", "sex")

# Create an empty list to store the data_coverage data frames
data_coverage_list <- list()

# Write a loop function
for (covariate in additional_covariates) {
  # Create Survival object 
  surv_obj <- Surv(time = data_surv$tte, event = data_surv$status)
  # Fit survival model for each ethnicity, IMD and covariates
  fit <- survfit(as.formula(paste("surv_obj ~ ethnicity + imd_Q5 +", covariate)), data = data_surv)
  # Plot the survival curves 
  ggsurv_obj <- ggsurvplot(fit, data = data_surv, risk.table = TRUE)
  # Extract the plot from ggsurv_obj, then save the plot
  ggsurv_obj$plot <- ggsurv_obj$plot + theme(legend.position = "none")
  ggsave(file.path(outdir, paste0("km_plot_all_groups_", covariate, ".png")), ggsurv_obj$plot)
  
  # Create the survival table
  surv_table <- ggsurv_obj$data.survplot
  # Filter to keep time for 12 weeks and 26 weeks
  data_coverage <- surv_table %>%
    filter(time %in% c(12*7, 26*7)) %>%
    mutate(coverage = 1 - surv) %>%
    select(strata, time, coverage, n.event, n.censor, std.err)
  
  # Add the data_coverage data frame to the list
  data_coverage_list[[covariate]] <- data_coverage
}

# Bind all the data_coverage data frames together
data_coverage_all <- bind_rows(data_coverage_list, .id = "covariate")

# Write the data to a CSV file
write_csv(
  data_coverage_all,
  file.path(outdir, glue::glue("vaccine_coverage_all_midpoint{threshold}.csv"))
)

##########
##########

# Create 5 plots: one for each 

# List of unique ethnicities
ethnicities <- unique(data_surv$ethnicity)
# Create a plot for each ethnicity
for (ethnicity in ethnicities) {
  data_subset <- data_surv[data_surv$ethnicity == ethnicity, ]
  surv_obj <- Surv(time = data_subset$tte, event = data_subset$status)
  fit <- survfit(surv_obj ~ imd_Q5, data = data_subset)
  plot <- ggsurvplot(fit, data = data_subset, risk.table = TRUE,
                     title = paste("Kaplan-Meier plot for", ethnicity),
                     xlab = "Time (days)", ylab = "Survival probability",
                     legend.title = "IMD Quintile",
                     palette = "jco",
                     theme = theme_classic2())
  # Save the plot to a file
  ggsave(file.path(outdir, paste0("km_plot_", gsub(" ", "_", ethnicity), ".png")), plot$plot)
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

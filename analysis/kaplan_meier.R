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

# Extract the plot from ggsurv_ob, then save the plot
ggsurv_obj$plot <- ggsurv_obj$plot + theme(legend.position = "none")
ggsave(file.path(outdir, "km_plot_all_groups.png"), ggsurv_obj$plot)

# Create a new dataset for vaccine coverage

# Create the survival table
surv_table <- ggsurv_obj$data.survplot %>%
  # Group by the specified covariates
  group_by(ethnicity, imd_Q5) %>%
  # Calculate the cumulative number of events and censors
  mutate(
    cum_n.event = cumsum(n.event),
    cum_n.censor = cumsum(n.censor)
  ) %>%
  ungroup()

# Filter the survival table keep time for 12 weeks and 26 weeks
data_coverage <- surv_table %>%
  filter(time %in% c(12*7, 26*7)) %>%
  mutate(
    coverage = 1 - surv,
    coverage.lower = 1 - upper,
    coverage.upper = 1 - lower
  ) %>%
  select(
    ethnicity, imd_Q5, 
    time, coverage, cum_n.event, cum_n.censor, std.err, 
    coverage.lower, coverage.upper 
  )

# Create a list of covariates
additional_covariates <- c("region", "jcvi_group", "sex")

# Create an empty list to store the data_coverage data frames
data_coverage_list <- list()

# Loop function
for (covariate in additional_covariates) {
  # Create Survival object inside the loop
  surv_obj <- Surv(time = data_surv$tte, event = data_surv$status)
  
  # Fit survival model for each ethnicity, IMD, and covariates
  fit <- survfit(as.formula(paste("surv_obj ~ ethnicity + imd_Q5 +", covariate)), data = data_surv)
  
  # Plot the survival curves 
  ggsurv_obj <- ggsurvplot(fit, data = data_surv, risk.table = TRUE)
  
  # Create the survival table
  surv_table <- ggsurv_obj$data.survplot %>%
    # Group by the specified covariates
    group_by(ethnicity, imd_Q5, !!sym(covariate)) %>%
    # Calculate the cumulative number of events and censors within groups
    mutate(
      cum_n.event = cumsum(n.event),
      cum_n.censor = cumsum(n.censor)
    ) %>%
    ungroup() 
  
  # Filter to keep time for 12 weeks and 26 weeks
  data_coverage_list[[covariate]] <- surv_table %>%
    filter(time %in% c(12*7, 26*7)) %>%
    mutate(
      coverage = 1 - surv,
      coverage.lower = 1 - upper,
      coverage.upper = 1 - lower
    ) %>%
    # Save the covariate as a column in the dataset
    mutate(covariate = covariate) %>%
    # Rename the covariate column to "level"
    rename(level = !!sym(covariate)) %>%
    select(
      ethnicity, imd_Q5, covariate, level, 
      time, coverage, cum_n.event, cum_n.censor, std.err,
      coverage.lower, coverage.upper
    )
}

# Bind all the data_coverage data frames together
data_coverage_all <- bind_rows(data_coverage_list, .id = "covariate")
# also bind to the previous one so you can save them all in one file to make
# it easier for the output checker
data_coverage_all <- bind_rows(
  data_coverage, # results for just ethnicity and imd
  data_coverage_all # results with additional variables
) %>% as_tibble()

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

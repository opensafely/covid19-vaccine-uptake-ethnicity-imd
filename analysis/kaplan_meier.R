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
    # combine the ethnicity and IMD categories - we will calculate kaplan-meier
    # estimates in each of these categories - you can do this by adding this 
    # variable to the right-hand side of the formula in survfit
    ethnicity_imd = paste0(ethnicity, ", ", imd_Q5),
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
fit_ethnicity_imd <- survfit(surv_obj ~ ethnicity_imd, data = data_surv)

# Plot the survival curves for each ethnicity and IMD subgroup
ggsurv_obj <- ggsurvplot(fit_ethnicity_imd, data = data_surv, risk.table = TRUE)

# Extract the plot from ggsurv_obj and get rid of the legend using `+ ggplot2::theme()`, then save the plot
ggsurv_obj$plot <- ggsurv_obj$plot + theme(legend.position = "none")
ggsave(file.path(outdir, "km_plot_all_groups.png"), ggsurv_obj$plot)

# Create a new dataset for vaccine coverage

# Create the survival table
surv_table <- ggsurv_obj$data.survtable 
# Calculate the survival probability and coverage
surv_table <- surv_table %>% # the object surv_table doesn't exist?
  mutate(surv = (n.risk - n.event) / n.risk,
         coverage = 1 - surv)
# Filter the survival table to only keep time = 12 weeks or 26 weeks
data_coverage <- surv_table %>%
  filter(time %in% c(12*7, 26*7)) %>%
  select(strata, time, coverage)
# Write the data to a CSV file
write_csv(
  data_coverage,
  file.path(outdir, glue::glue("vaccine_coverage_midpoint{threshold}.csv"))
)

# Create 5 plots: one for each 

# List of unique ethnicities
ethnicities <- unique(sapply(data_surv$ethnicity_imd, function(x) strsplit(x, ", ")[[1]][1]))

# Create a plot for each ethnicity
for (ethnicity in ethnicities) {
  # Subset the data for the current ethnicity
  data_subset <- data_surv[grep(ethnicity, data_surv$ethnicity_imd), ]
  # Create a Surv object for kaplan-meier analysis
  surv_obj <- Surv(time = data_subset$tte, event = data_subset$status)
  # Fit a survival curve
  fit <- survfit(surv_obj ~ ethnicity_imd, data = data_subset)
  # Plot the survival curve and save it to a variable
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

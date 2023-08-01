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

# # Fit the Kaplan-Meier survival model
# fit <- survfit(surv_obj ~ 1)
# 
# # Plot the survival curve
# ggsurvplot(fit, data = data_surv, risk.table = TRUE)

#########
#########


# Fit the Kaplan-Meier survival model for each ethnicity and IMD subgroup
fit_ethnicity_imd <- survfit(surv_obj ~ ethnicity_imd, data = data_surv)

# Plot the survival curves for each ethnicity and IMD subgroup
ggsurv_obj <- ggsurvplot(fit_ethnicity_imd, data = data_surv, risk.table = TRUE)

# TODO
# 1) extract the plot from ggsurv_obj and get rid of the legend using 
#    `+ ggplot2::theme()`, then save the plot
# 2) create a new dataset by filtering the data to only keep time = 12*7 or 26*7, 
#    and use 1-surv to calculate the vaccine coverage. only keep the columns 
#    that you wan to release from opensafely, and save this dataset as a csv file.
# 3) update your code in the for loop below to create 5 plots: one for each 
#    ethnicity where the IMD_Q5 groups have different colour lines, save these 
#    plots.

##########
##########

# List of unique ethnicity and IMD subgroups
ethnicity_imd_groups <- unique(data_surv$ethnicity_imd)

# Create a plot for each ethnicity and IMD subgroup
for (group in ethnicity_imd_groups) {
  # Subset the data for the current group
  data_subset <- data_surv[data_surv$ethnicity_imd == group, ]
  
  # Create a Surv object for kaplan-meier analysis
  surv_obj <- Surv(time = data_subset$tte, event = data_subset$status)
  
  # Fit a survival curve
  fit <- survfit(surv_obj ~ 1, data = data_subset)
  
  # Plot the survival curve and save it to a variable
  plot <- ggsurvplot(fit, data = data_subset, risk.table = TRUE,
                     title = paste("Kaplan-Meier plot for", group),
                     xlab = "Time (days)", ylab = "Survival probability",
                     legend.title = "Status",
                     palette = c("#E7B800", "#2E9FDF"),
                     theme = theme_classic2())
  
  # Save the plot to a file
  ggsave(paste0("km_plot_", gsub(" ", "_", gsub(",", "", group)), ".png"), plot$plot)
}



# Check unique values of the 'status' variable
unique(data_surv$status)




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

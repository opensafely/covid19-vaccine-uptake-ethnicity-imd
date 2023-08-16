
library(tidyverse)
library(here)
library(viridis)

# Source metadata and functions
source(here("analysis", "design.R"))

outdir <- here("release20230809")

# Read the data 
data_vax_counts <- readr::read_csv(
  file.path(outdir, glue("data_vax_counts_midpoint{threshold}.csv"))
) %>%
  mutate(across(elig_date_rank, ~as.factor(.x)))

# Plot using geom_line()
data_vax_counts %>%
  ggplot(aes(x = covid_vax_disease_1_time, y = n, color = elig_date_rank)) +
  geom_line() +
  facet_wrap(~jcvi_group, scales = "free_y", nrow = 4) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(nrow = 1)) +
  theme_minimal() +
  theme(legend.position = c(0.65, 0.1)) +
  labs(x = "Days between eligibility and first vaccination",
       y = "Number of patients",
       color = "Eligibility Date")
# Save
ggsave(
  file.path(outdir, "vax_dates_line.png")#,
  # width, height, units = "cm"
  )



#############
#############



# Read the data 
data_bar_plots <- readr::read_csv(
  file.path(outdir, glue("combined_data_bar_plot_midpoint{threshold}.csv"))
)
data_bar_plots


# Function for bar plot
create_bar_plot <- function(data, path) {
  variable_name <- unique(data$variable)
  levels <- unique(data$level)
  # fill_pal <- Polychrome::light.colors(n=24) # Polychrome not installed in OpenSAFELY
  fill_pal <- c(
    "#FD3216", "#00FE35", "#6A76FC", "#FED4C4", "#FE00CE", "#0DF9FF", "#F6F926",
    "#FF9616", "#479B55", "#EEA6FB", "#DC587D", "#D626FF", "#6E899C", "#00B5F7",
    "#B68E00", "#C9FBE5", "#FF0092", "#22FFA7", "#E3EE9E", "#86CE00", "#BC7196", 
    "#7E7DCD", "#FC6955", "#E48F72"
  )
  if (length(levels) == 2) {
    # otherwise it picks bright red and green which doens't look nice...
    fill_pal <- fill_pal[c(24,14)]
  } else {
    # you could pick a different colour scheme from Polychrome or another package if you prefer
    fill_pal <- fill_pal[-c(1:2)] # this is personal preference, I just don't like the red and green
    fill_pal <- fill_pal[1:length(levels)]
  }
  names(fill_pal) <- levels
  p <- ggplot(data, aes(x = level, y = percent, fill = level)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = fill_pal) +
    facet_grid(rows = vars(ethnicity), cols = vars(imd_Q5)) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) +
    labs(x = NULL,
         y = "Percentage (%)",
         fill = variable_name)
  ggsave(
    filename = file.path(path, paste0("strat_dist_", variable_name, ".png")),
    plot = p
  )
  return(p)
}

list_of_data <- split(data_bar_plots, data_bar_plots$variable)



plots <- lapply(
  seq_along(list_of_data), 
  function(x) create_bar_plot(list_of_data[[x]], path = outdir)
)



#############
#############

# Read the data 
vaccine_coverage <- readr::read_csv(
  file.path(outdir, glue("vaccine_coverage_all_midpoint{threshold}.csv"))
)

vaccine_coverage

# Generate the plot for vaccine coverage across different ethnicities and IMD subgroups at 12 and 26 weeks.
ggplot(vaccine_coverage, aes(x = ethnicity, y = coverage, fill = imd_Q5, group = interaction(ethnicity, imd_Q5))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(~time) + 
  labs(
    title = "Vaccine Coverage Across Ethnicities and IMD subgroups",
    x = "Ethnicity",
    y = "Coverage Percentage",
    fill = "IMD Q5 Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#########
#########


 
# List of unique ethnicities
ethnicities <- unique(vaccine_coverage$ethnicity)
# Create a plot for each ethnicity
for (ethnicity in ethnicities) {
   data_subset <- vaccine_coverage[vaccine_coverage$ethnicity == ethnicity, ]
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























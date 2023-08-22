# install.packages("RColorBrewer")


library(survival)
library(survminer)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(here)
library(viridis)
library(dplyr)
library(tidyr)
library(grid)

# Source metadata and functions
source(here("analysis", "design.R"))

outdir <- here("release20230809")

#############

# Define the jcvi_groups
jcvi_groups <- 
  tribble(
    ~group, ~definition,
    "01", "longres And age > 65",
    "02", "age >=80",
    "03", "age >=75",
    "04a", "age >=70",
    "04b", "cev_group AND age >=18", 
    "05", "age >=65",
    "06", "atrisk_group AND age >=18", 
    "07", "age >=60",
    "08", "age >=55",
    "09", "age >=50",
    "10", "age >=40",
    "11", "age >=30",
    "12", "age >=18",
    "99", "DEFAULT",
  )

# Read the data 
data_vax_counts <- readr::read_csv(
  file.path(outdir, glue("data_vax_counts_midpoint{threshold}.csv"))
) %>%
  mutate(across(elig_date_rank, ~as.factor(.x)))

# Join with jcvi_groups to get the definition for each group
data_vax_counts <- left_join(data_vax_counts, jcvi_groups, by = c("jcvi_group" = "group"))

# Plot using geom_line()
p <- data_vax_counts %>%
  ggplot(aes(x = covid_vax_disease_1_time, y = n, color = elig_date_rank)) +
  geom_line(size=1) + 
  # Use the definition variable for the facet labels
  facet_wrap(~definition, scales = "free_y", nrow = 4) +
  scale_color_viridis_d(name="Eligibility Date Rank") +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.position = "bottom",
    legend.title = element_text(face="bold"),
    plot.title = element_text(face="bold", hjust=0.5),
    axis.title.x = element_text(face="bold"),
    axis.title.y = element_text(face="bold")
  ) +
  labs(
    title="COVID-19 Vaccination Uptake Over Time",
    x = "Days Post-Eligibility",
    y = "Number of Patients Vaccinated"
  )

# Save the plot 
ggsave(
  file.path(outdir, "vax_dates_line.png"),
  width = 10, 
  height = 8, 
  units = "in", 
  dpi = 300 
)



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

# Function for bar plot
create_bar_plot <- function(data, path) {
  variable_name <- unique(data$variable)
  levels <- unique(data$level)
  
  # color palette 
  fill_pal <- scales::hue_pal()(length(levels))
  names(fill_pal) <- levels
  
  p <- ggplot(data, aes(x = level, y = percent, fill = level)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) + 
    scale_fill_manual(values = fill_pal) +
    facet_grid(rows = vars(ethnicity), cols = vars(imd_Q5)) +
  
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      #legend.position = "bottom",
      plot.title = element_text(face="bold", hjust=0.5),
      axis.title.x = element_text(face="bold"),
      axis.title.y = element_text(face="bold"),
      axis.text.x = element_blank() # Remove x-axis labels
    ) +
    labs(
      title = paste("Distribution of", variable_name, "by Ethnicity and IMD Quintile"),
      x = NULL,
      y = "Percentage (%)"
    )
  
  # Save the plot 
  ggsave(
    filename = file.path(path, paste0("strat_dist_", variable_name, ".png")),
    plot = p,
    width = 10, 
    height = 8, 
    units = "in", 
    dpi = 300 
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

# Function for dot plot
create_dot_plot <- function(data, path) {
  variable_name <- unique(data$variable)
  levels <- unique(data$level)
  
  # Use a more muted color palette
  color_pal <- colorRampPalette(c("lightblue", "darkblue"))(length(levels))
  names(color_pal) <- levels
  
  p <- ggplot(data, aes(x = level, y = percent, color = level)) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.7) +  # Add jitter to avoid overplotting
    facet_grid(rows = vars(ethnicity), cols = vars(imd_Q5)) +
    scale_color_manual(values = color_pal, name = variable_name) +
    theme_bw() +  # Use a cleaner theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14)
    ) +
    labs(
      x = NULL,
      y = "Percentage (%)",
      title = "COVID-19 Vaccination Uptake by Demographic Groups",
      subtitle = paste("Stratified by", variable_name)
    )
  
  ggsave(
    filename = file.path(path, paste0("dot_plot_", variable_name, ".png")),
    plot = p,
    width = 10, 
    height = 6,
    units = "in", 
    dpi = 300
  )
  return(p)
}

list_of_data <- split(data_bar_plots, data_bar_plots$variable)

# Create dot plots
dot_plots <- lapply(
  seq_along(list_of_data), 
  function(x) create_dot_plot(list_of_data[[x]], path = outdir)
)
#############
#############

#############
#############

# Read the data 
data_surv <- readr::read_csv(
  file.path(outdir, glue("vaccine_coverage_all_midpoint{threshold}.csv"))
)

# List of unique ethnicities
ethnicities <- unique(data_surv$ethnicity)

# Create a plot for each ethnicity
for (ethnicity in ethnicities) {
  data_subset <- subset(data_surv, ethnicity == ethnicity)
  
  # Create a ggplot for vaccine coverage across different IMD subgroups
  plot <- ggplot2::ggplot(data_subset, ggplot2::aes(x = time, y = coverage, color = as.factor(imd_Q5))) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::labs(title = paste("Vaccine Coverage for", ethnicity),
                  x = "Time (days)",
                  y = "Coverage",
                  color = "IMD Quintile") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title.x = ggplot2::element_text(face = "bold", size = 14),
      axis.title.y = ggplot2::element_text(face = "bold", size = 14),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "bottom",
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_brewer(palette = "Set1")
  
  
  # Save the plot 
  ggplot2::ggsave(file.path(outdir, paste0("vaccine_coverage_", gsub(" ", "_", ethnicity), ".png")), plot, width = 10, height = 6, dpi = 300)
}
###########
###########



data_surv %>%
  filter(covariate == 'jcvi_group', level %in% c("10","11","12")) %>%   #grouping for younger population
  ggplot(
    aes(
      x = time, y = coverage,
      colour = imd_Q5
    )
  ) +
  geom_errorbar(
    aes(ymin = coverage.lower, ymax = coverage.upper),
    width = 10
    ) +
  geom_point(size = 1) +
  geom_line() +
  facet_grid(col = vars(ethnicity)) + #rows = vars(ethnicity)
  scale_x_continuous(
    breaks = c(84,182)
  ) +
  scale_colour_viridis_d() + 
  labs(
    x = "Days since eligible",
    title = "Vaccine Coverage in London by Ethnicity and IMD Quintile"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

# Save the plot 
ggsave(
  file.path(outdir, "vax_coverage_title.png"),
  width = 10, 
  height = 8, 
  units = "in", 
  dpi = 300 
)
###########
###########

data_surv %>%
  filter(covariate == 'jcvi_group', level %in% c("01","02","03","04","05")) %>%
  group_by(time, imd_Q5, ethnicity) %>%
  summarise(
    mean_coverage = mean(coverage, na.rm = TRUE),
    coverage_lower = mean(coverage.lower, na.rm = TRUE),
    coverage_upper = mean(coverage.upper, na.rm = TRUE),
    .groups = 'drop' # This line takes care of the warning about grouped output
  ) %>%
  ggplot(
    aes(
      x = time, y = mean_coverage,
      colour = imd_Q5
    )
  ) +
  geom_errorbar(
    aes(ymin = coverage_lower, ymax = coverage_upper),
    width = 10
  ) +
  geom_point(size = 1) +
  geom_line() +
  facet_grid(col = vars(ethnicity)) +
  scale_x_continuous(
    breaks = c(84,182)
  ) +
  scale_colour_viridis_d() + 
  labs(
    x = "Days since eligible",
    title = "Vaccine Coverage in Adults >65 by Ethnicity and IMD Quintile"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(outdir, "vax_coverage_Age_65.png"),
  width = 10, 
  height = 8, 
  units = "in", 
  dpi = 300 
)
################
################

summary_by_ethnicity <- data_surv %>%
  group_by(ethnicity) %>%
  summarise(
    mean_coverage = mean(coverage, na.rm = TRUE),
    sd_coverage = sd(coverage, na.rm = TRUE),
    min_coverage = min(coverage, na.rm = TRUE),
    max_coverage = max(coverage, na.rm = TRUE),
    mean_coverage_lower = mean(coverage.lower, na.rm = TRUE),
    mean_coverage_upper = mean(coverage.upper, na.rm = TRUE),
    .groups = 'drop' 
  )

summary_by_ethnicity









# TODO
# split by the covariates adding them as rows in the facets
# you might want to do this a few levels at a time for jcvi_group and region,
# otherwise the plot will be quite hard to read


# Create the plot
plot <- data_surv %>%
  filter(is.na(covariate)) %>%
  ggplot(
    aes(
      x = time, y = coverage,
      colour = imd_Q5
    )
  ) +
  geom_errorbar(
    aes(ymin = coverage.lower, ymax = coverage.upper),
    width = 10
  ) +
  geom_point(size = 1) +
  geom_line() +
  facet_grid(rows = vars(jcvi_group), cols = vars(ethnicity)) + # Faceting by jcvi_group
  scale_x_continuous(
    breaks = c(84,182)
  ) +
  scale_colour_viridis_d() + 
  labs(
    x = "Days since eligible"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )
 
print(plot)



############
############

#flowchart summary table

flowchart <- readr::read_csv(
  file.path(outdir, glue("flowchart_midpoint{threshold}.csv"))
)


sum_table <- flowchart[, c("criteria", "n", "n_exclude", "pct_exclude")]

# Convert the percentage 
sum_table$pct_exclude <- sum_table$pct_exclude * 100

# Convert the data frame to a table grob
table_grob <- gridExtra::tableGrob(sum_table, rows = NULL) 

# Save the table as a PNG image
png(filename = paste0(outdir, "/sum_table.png"), width = 800, height = 600)
grid.draw(table_grob)
dev.off()




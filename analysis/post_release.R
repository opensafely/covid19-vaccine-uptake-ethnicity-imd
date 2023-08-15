
library(tidyverse)
library(here)
library(viridis)

# Source metadata and functions
source(here("analysis", "design.R"))

outdir <- here("release20230809")

# Save to .csv file for release
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

library(tidyverse)
library(here)

# Source metadata and functions
source(here("analysis", "design.R"))

outdir <- here("release20230809")

# Read the data 

data_surv <- readr::read_csv(
  file.path(outdir, glue("vaccine_coverage_all_midpoint{threshold}.csv"))
)

data_bar_plots <- readr::read_csv(
  file.path(outdir, glue("combined_data_bar_plot_midpoint{threshold}.csv"))
)

data_surv %>%
  filter(covariate == "sex") %>%
  arrange(ethnicity, imd_Q5, covariate, time) %>%
  group_by(ethnicity, imd_Q5, covariate, time) %>%
  mutate(p.cum_n.event = 100*cum_n.event/sum(cum_n.event)) %>%
  left_join(
    data_bar_plots,
    by = c("ethnicity", "imd_Q5", "covariate" = "variable", "level")
    ) %>%
  select(
    ethnicity, imd_Q5, covariate, level, time,
    `Percent of vaccines` = p.cum_n.event, 
    `Percent of people` = percent
    ) %>%
  pivot_longer(
    cols = starts_with("Percent")
  ) %>%
  filter(time == 182) %>%
  # filter(time == 84) %>%
  ggplot(aes(x = level, y = value, fill = name)) +
  geom_hline(yintercept = 50, colour = "grey", linetype = "dashed") +
  # geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  geom_point(shape = 21, alpha = 0.75, colour=alpha("black", 1)) +
  facet_grid(cols = vars(ethnicity), rows = vars(imd_Q5)) +
  scale_fill_discrete(name = NULL) + 
  labs(x = "Sex", y = "Percent") +
  theme_bw() +
  theme(legend.position = "bottom") #+
  # coord_cartesian(ylim = c(44,56))
  
ggsave(
  file.path(outdir, paste0("charlotte_suggestion.png")), 
  width = 10, height = 7, dpi = 300
  )

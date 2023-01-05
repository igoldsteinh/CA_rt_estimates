# visualize results
library(tidyverse)
library(tidybayes)
library(gridExtra)
library(ggplot2)
library(scales)
library(cowplot)

data <- read_csv(here::here("data", "cases_hospitalizations_by_county.csv"))

county_region_key <- read_csv(here::here("data", "county_region_key.csv"))

rt_results <- read_csv(here::here("full_county_rt_estimates.csv")) 

rt_results <- rt_results %>% 
              left_join(county_region_key, by = "county")


rt_results$region[rt_results$county == "California"] <- "A"
# all credit to Damon Bayer for plot functions 
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))

make_rt_plot <- function(target_place_name) {
  rt_results %>%
    filter(county == target_place_name) %>%
    filter(!is.na(date)) %>%
    ggplot(aes(date, rt, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_y_continuous("Rt", label = comma) +
    scale_x_date(name = "Date") +
    ggtitle(str_c(target_place_name %>% str_to_title(), ifelse(target_place_name != "California", "County", ""), sep = " "),
            subtitle = str_c("Estimated", max(rt_results$date), sep = " ")) +
    my_theme
}

ggsave2(filename = here::here("figures", "rt_plots.pdf"),
        plot = rt_results %>%
          distinct(region, county) %>%
          arrange(region, county) %>%
          pull(county) %>%
          map(make_rt_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)

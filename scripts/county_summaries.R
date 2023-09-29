# calculates stan summary diagnostics and final output file for all CA counties
# and statewide estimate as well
library(tidyverse)
library(tidybayes)
library(stringr)
library(gridExtra)
library(ggplot2)
library(scales)
library(cowplot)

source(here::here("src", "rt_functions.R"))


# read in data and results ---------------------------------------------------------

# county_names <- c("Alameda")
# file_names <- c("sender_model.rds")
county_id_key <- read_csv("data/county_id_key.csv")
county_names <- county_id_key$county
county_id <- county_id_key$id
# county_names <- c( "California")
# county_id = 0


# then we have to reduce down the desired chains --------------------------
full_stan_diag <- map2(county_names, county_id, ~read_csv(here::here("results", 
                                                           "standiags",
                                                           str_c("id=", 
                                                                 .y, 
                                                                 "_", 
                                                                 .x, 
                                                                 "_standiags", 
                                                                 ".csv", 
                                                                 "")))) %>%
                  bind_rows()
print(full_stan_diag)
write_csv(full_stan_diag, here::here("results", "full_stan_diag.csv"))
# create final rt frame ---------------------------------------------------
full_county_rt <- map2(county_names, county_id, ~read_csv(here::here("results",
                                                         "rt_credible_intervals",
                                                         str_c("id=", 
                                                               .y, 
                                                               "_", 
                                                               .x, 
                                                               "_rt_intervals", 
                                                               ".csv", "")))) %>%
             bind_rows() %>%
             arrange(county) %>%
             dplyr::select(date, county, rt, .lower, .upper, .width, .point)

county_rt <- full_county_rt %>%
             filter(.width == 0.95) 
print("hello")
write_csv(full_county_rt, "full_county_rt_estimates.csv")
write_csv(county_rt, "CDPH_county_rt_estimates.csv")


data <- read_csv(here::here("data", "cases_hospitalizations_by_county.csv"))

county_region_key <- read_csv(here::here("data", "county_region_key.csv"))

rt_results <- full_county_rt

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
            subtitle = str_c("Estimated", max(rt_results$date, na.rm = TRUE), sep = " ")) +
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


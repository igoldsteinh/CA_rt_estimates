# calculates stan summary diagnostics and final output file for all CA counties
# and statewide estimate as well
library(tidyverse)
library(tidybayes)
library(stringr)
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




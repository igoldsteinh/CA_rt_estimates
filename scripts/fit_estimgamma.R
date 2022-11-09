#estimgamma fit to all california counties and statewide
library(tidyverse)
library(rstan)
library(tidybayes)
library(patchwork)
library(EpiEstim)
#library(epidemia)
#library(rstanarm)
library(lubridate)
#library(coda)
library(fs)
source(here::here("src", "rt_functions.R"))
dir_create(path("results", "posteriors"))
dir_create(path("results", "rt_credible_intervals"))


# command args for array job ----------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
indic <- as.integer(args[1])


# code start ---------------------------------------------------------------


set.seed(225)
options(mc.cores = parallelly::availableCores())
rstan_options(auto_write = TRUE)


ca_data <- read_csv("data/cases_hospitalizations_by_county.csv")

county_id_key <- read_csv("data/county_id_key.csv")

county_name <- county_id_key %>% filter(id == indic) %>% pull(county)
#create weekly data, currently using data from the week of 08/02/2020 - 01/09/2022
county_data <- ca_data %>% filter(county == county_name) %>% 
  rename(total_cases = est_cases, 
         total_tests = est_tests)

# read in priors for overdispersion
overdisp_filename <- paste0("overdisp_priors_countyid", indic, ".csv")
overdisp_prior <- read_csv(here::here("data", overdisp_filename))

# calculate quantile for tests
test_quantile <- quantile(county_data$total_tests)

# fit model

# choose starting points

# first choose rt starting points using epiestim
logrt_start <- get_logrtstart(county_data)

#next choose incidence starting points
incid_start <- 1/0.2 * county_data$total_cases

init_func <- function() list(log_incid_rate_raw = 0,
                             log_rt0_raw = 0,
                             rho = 0.2/test_quantile[2],
                             kappa = county_kappa$mean[1],
                             seed_incid_one_raw =1,
                             incid = incid_start,
                             log_rt = logrt_start)
# gen params are the generation time parameters, first one is the latent period rate
# second is the infectious period rate (scaled by 7 becuase we're using weeks)
# delay_params are for a gamma, right now I'm having it just be the latent period

# rho is the "case detection prior" rho * tests * cases is the mean of case observation model
# if you're fitting CA data, the rho prior should work well as is

# there are other priors but their defaults should be fine
county_posterior <- fit_estimgamma_model(county_data,
                                          gen_params = c(log(7.872346) + log(1/7), 
                                                         0.642713),
                                          delay_params = c(4.05, 7*0.74),
                                          prev_vals = 4,
                                          log_rho_mean = log(0.066/test_quantile[2]),
                                          log_rho_sd = 0.3,
                                          kappa_mean = county_kappa$mean[1],
                                          kappa_sd = county_kappa$sd[2],
                                         log_r0_mean = log(1.5),
                                         log_r0_sd = 0.75,
                                          iterations = 100,
                                          init_func = init_func,
                                         gen_dist = "log-normal",
                                         seed = 225,
                                         thin = 3)


rt_posterior_summary  <- summarise_realdata_rt_estimgamma(stan_posterior = county_posterior,
                                                          weekly_data = county_data, 
                                                          start_date = 1,
                                                          include_chains = c(1,2,3)) %>%
                        mutate(county = county_name) 

write_rds(county_posterior, here::here("results", 
                                      "posteriors",
                                      str_c("id=", indic, "_", county_name, "_estimgamma_posterior", ".rds", "")))

write_csv(rt_posterior_summary, here::here("results",
                                           "rt_credible_intervals",
                                           str_c("id=", indic, "_", county_name, "_rt_intervals", ".csv", "")))
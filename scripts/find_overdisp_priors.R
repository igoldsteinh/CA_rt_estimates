# find parameters for overdispersion priors
# based on posteriors from NB spline
library(tidyverse)
library(rstan)
library(brms)
options(mc.cores = parallelly::availableCores())
rstan_options(auto_write = TRUE)

source("src/rt_functions.R")

# command args for array job ----------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
indic <- as.integer(args[1])

county_id_key <- read_csv("data/county_id_key.csv")
county_id <- county_id_key %>% filter(id == indic) %>% pull(county)

county_data <- read_csv("data/cases_hospitalizations_by_county.csv") %>%
        filter(county == county_id) %>% 
  rename(total_cases = est_cases, 
         total_tests = est_tests)

# run spline for kappa priors
county_spline <- run_nb_spline(county_data)

# calculate kappa priors
county_kappa <- choose_kappa_params(county_spline)

# save the results --------------------------------------------------------
labels <- c("cases")

priors <- t(county_kappa$par) %>%
          cbind(labels)

priors <- data.frame(priors)

colnames(priors) <- c("mean", "sd", "labels")
rownames(priors) <- NULL

file_name <- paste("overdisp_priors_countyid", indic, ".csv", sep = "")
write_csv(priors, paste("data/", file_name, sep = ""))

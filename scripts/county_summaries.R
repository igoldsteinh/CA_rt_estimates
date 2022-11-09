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
# county_names <- c("Alameda", "California")

county_posteriors <- map(county_names, ~read_rds(here::here("results", "posteriors", paste0(.x, "_estimgamma_posterior.rds"))))
# first task just look at the traces --------------------------------------
num_counties <- length(county_names)

for (i in 1:num_counties) {
  
  trace <- rstan::traceplot(county_posteriors[[i]], pars = "lp__") +
    ggtitle(county_names[i])
  
  ggsave(here::here("figures", str_c(county_names[i], "_trace", ".pdf", sep = "")), 
         plot = trace, 
         width = 5, 
         height = 5)
  

  
}



# then we have to reduce down the desired chains --------------------------
good_chains = c(1,2,3)
draws <- map(county_posteriors, ~.x %>% tidy_draws())

standiags <- map(draws, calc_stan_diags, include_chains = good_chains)


rhat_frame <- map(standiags, pluck, 2) %>%
              map(data.frame) %>%
              map2(county_names, ~.x %>%
                     mutate(county = .y,
                            rownumber = row_number()) %>%
                     rename(rhat = 1)) %>%
              bind_rows()

bulkess_frame <-  map(standiags, pluck, 3) %>%
  map(data.frame) %>%
  map2(county_names, ~.x %>% mutate(county = .y,
                                    rownumber = row_number()) %>%
                       rename(bulkess = 1)) %>%
  bind_rows()

tailess_frame <- map(standiags, pluck, 4) %>%
  map(data.frame) %>%
  map2(county_names, ~.x %>% mutate(county = .y,
                                    rownumber = row_number()) %>%
                      rename(tailess = 1)) %>%
  bind_rows()

rt_ess <- map(standiags, pluck, 1) %>%
          map(data.frame) %>%
          map(~.x %>%
                filter(str_detect(var_names, "log_rt")) %>%
                ungroup() %>%
                summarise(min_rt_essbulk = min(essbulk),
                          min_rt_esstail = min(esstail))) %>%
          map2(county_names, ~.x %>% mutate(county = .y,
                                            rownumber = 1)) %>%
          bind_rows()

# print(names(rt_ess))
# print(names(bulkess_frame))
# print(names(rhat_frame))
# print(names(tailess_frame))
standiag_frame <- rhat_frame %>%
                  left_join(bulkess_frame, by = c("county" = "county",
                                                  "rownumber"="rownumber")) %>%
                  left_join(tailess_frame, by= c("county" = "county",
                                                 "rownumber" = "rownumber") ) %>%
                  left_join(rt_ess, by= c("county" = "county",
                                          "rownumber" = "rownumber"))

# add in bfmi
bfmi <- map(county_posteriors, get_bfmi)
bfmi_range <- map(bfmi, ~t(range(.x))) %>%
              map(data.frame) %>%
              bind_rows() %>%
              rename(min_bfmi = 1, 
                     max_bfmi = 2) %>%
              cbind(county_names)
# print(bfmi)
standiag_frame <- standiag_frame %>%
                  left_join(bfmi_range, by = c("county" = "county_names"))


write_csv(standiag_frame, here::here("results", "ca_county_standiags.csv"))



# create final rt frame ---------------------------------------------------
full_county_rt <- map(county_names, ~read_csv(here::here("results", 
                                                            "rt_credible_intervals", 
                                                            paste0(.x, "_rt_intervals.csv")))) %>%
             bind_rows() %>%
             arrange(county)

county_rt <- full_county_rt %>%
             filter(.width == 0.95)

write_csv(full_county_rt, "full_county_rt_estimates.csv")
write_csv(county_rt, "CDPH_county_rt_estimates.csv")




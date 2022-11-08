# damon functions
library(tidyverse)
library(EpiEstim)
library(tidybayes)
library(truncnorm)
library(rstan)
library(lubridate)
library(brms)
library(epidemia)
library(rstanarm)
library(sdprisk)

# Create automated process for choosing kappa parameter from a spl --------

run_nb_spline <- function(data, 
                          seed = 17,
                          iter = 4000,
                          warmup = 1000,
                          thin = 10, 
                          refresh = 0,
                          adapt_delta = 0.99) {
  spline_model <- brm(bf(total_cases ~ s(time)),
                      data = data, family = negbinomial(), cores = 4, seed = seed,
                      iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                      control = list(adapt_delta = adapt_delta))
  
  return(spline_model)
}

compare_kappa_quantiles <- function(candidate_params, true_quantiles) {
  candidate_quantiles <- quantile(rtruncnorm(10000, 
                                             a = 0, 
                                             mean = candidate_params[1], 
                                             sd = candidate_params[2]),
                                  c(0.025, 0.975))
  loss <- (true_quantiles[1] - candidate_quantiles[1])^2 + (true_quantiles[2] - candidate_quantiles[2])^2
  return(loss)
}


choose_kappa_params <- function(spline_posterior) {
  posterior_pars <- summary(spline_posterior)
  start_mean <- posterior_pars[["spec_pars"]][[1]]
  start_sd <- posterior_pars[["spec_pars"]][[2]]
  true_lb <- posterior_pars[["spec_pars"]][[3]]
  true_ub <- posterior_pars[["spec_pars"]][[4]]
  
  start_params <- c(start_mean, start_sd)
  true_quantiles <- c(true_lb, true_ub)
  
  optim_params <- optim(par = start_params,
                        fn = compare_kappa_quantiles, 
                        true_quantiles = true_quantiles)
}

# fit exp_seed model ------------------------------------------------------
# this is the current version of the model
# assuming generation time is a hypo-expo distribution
# assuming delay time is a gamma distribution, admitting the zero case
fit_estimgamma_model <- function(data,
                                 gen_params,
                                 delay_params,
                                 prev_vals,
                                 log_nu_mean = -2,
                                 log_nu_sd = 0.7,
                                 log_sigma_mean = -0.6,
                                 log_sigma_sd = 0.6,
                                 log_rho_mean,
                                 log_rho_sd,
                                 log_r0_mean = log(1),
                                 log_r0_sd = 0.75,
                                 kappa_mean,
                                 kappa_sd,
                                 init_func,
                                 iterations = 2000,
                                 thin = 2,
                                 adapt_delta = 0.99,
                                 treedepth = 12,
                                 seed = 45,
                                 chain = 4,
                                 gen_dist = "hypo-exp",
                                 delay_dist = "gamma") {
  
  data_length <- dim(data)[1]
  
  if (gen_dist == "hypo-exp") {
    gen_weights <- epidemia_hypoexp(data_length, gen_params)
    
  }
  
  
  
  if (gen_dist == "log-normal") {
    gen_weights <- epidemia_lognormal(data_length, gen_params)
  }
  
  if (gen_dist == "weibull") {
    gen_weights <- epidemia_weibull(data_length, gen_params)
  }
  
  
  if (delay_dist == "gamma") {
    delay_weights <- zero_epidemia_gamma(data_length, 
                                         delay_params[1], 
                                         delay_params[2])
  }
  
  model_object <- list(n = data_length, 
                       d = data_length,
                       w = gen_weights,
                       delay_weights = delay_weights,
                       obs = data$total_cases,
                       test = data$total_tests,
                       prev_vals = 4,
                       log_incid_rate_mean = log_nu_mean,
                       log_incid_rate_sd = log_nu_sd,
                       log_sigma_mu = log_sigma_mean,
                       log_sigma_sd = log_sigma_sd,
                       log_rho_mu = log_rho_mean,
                       log_rho_sd = log_rho_sd,
                       log_r0_mu = log_r0_mean,
                       log_r0_sd = log_r0_sd,
                       kappa_mu = kappa_mean,
                       kappa_sd = kappa_sd)
  
  
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = treedepth)
  
  model_fit <- stan(file = "src/rt_estim_gamma.stan",
                    data = model_object,
                    seed = seed,
                    iter = iterations,
                    thin = thin,
                    chain = chain,
                    init = init_func,
                    control = control_list)
  
  return(model_fit)
}

# create estimgamma rt posterior for sim data, no truth involved --------------------
summarise_rt_estimgamma <- function(stan_posterior,
                                             start_date,
                                             include_chains = c(1,2,3,4)){
  
  
  rt_posterior <- stan_posterior %>%
    spread_draws(log_rt[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i + 0 + start_date -1) %>%
    mutate(rt = exp(log_rt)) %>%
    dplyr::select(date, rt) %>%
    group_by(date) %>% 
    median_qi(.width = c(0.5, 0.8, 0.95)) 
  return(rt_posterior)
}


# rt_metrics --------------------------------------------------------------
# operating characteristics

rt_metrics<- function(data, value, upper, lower) {
  metric_one <- data %>%
    mutate(dev = abs({{ value }} - true_rt),
           CIW = abs({{ upper }} - {{ lower }}),
           envelope = true_rt >= {{ lower }} & true_rt <=  {{ upper }}) %>%
    ungroup() %>%
    filter(!is.na(dev)) %>%
    summarise(mean_dev = mean(dev),
              MCIW = mean(CIW),
              mean_env = mean(envelope))
  
  metrics_two <- data %>%
    mutate(prev_val = lag({{ value }}),
           prev_rt = lag(true_rt),
           sv = abs({{ value }} - prev_val),
           rt_sv = abs(true_rt - prev_rt)) %>%
    filter(!is.na(sv)) %>%
    ungroup() %>%
    summarise(MASV = mean(sv),
              true_MASV = mean(rt_sv))
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}


# use epiestim to choose initial conditions for rt ------------------------

get_logrtstart <- function(data,
                           window = 1, 
                           GI_mean = 11.5/7
) {
  
  window = window
  GI_mean = GI_mean
  GI_var = 2*(GI_mean/2)^2
  
  ts <- data$time
  ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
  te <- ts+(window-1)
  
  estimate_R(
    incid = data$total_cases,
    method = "uncertain_si",
    config = make_config(
      list(
        mean_si = GI_mean,
        min_mean_si = 1,
        max_mean_si = GI_mean + 1,
        std_mean_si = 1.5,
        std_std_si = 1.5,
        std_si = sqrt(GI_var),
        min_std_si = sqrt(GI_var)*.8,
        max_std_si = sqrt(GI_var)*1.2,
        n1 = 50,
        n2 = 100, 
        t_start=ts,
        t_end=te
      )
    )
  ) -> ee_outs
  
  ee_quantile <- ee_outs[["R"]] %>%
    dplyr::select(t_start, 
                  rt_mean = `Mean(R)`, 
                  rt_median = `Median(R)`,
                  rt_CI95l = `Quantile.0.025(R)`,
                  rt_CI95u = `Quantile.0.975(R)`) %>%
    mutate(time  = t_start) 
  
  
  log_ee_median <- log(ee_quantile %>% pull(rt_median))
  
  first_one <- log_ee_median[1]
  rt_start <- c(first_one, log_ee_median)
  
  return(rt_start)
}

# create rt posterior for real data, no truth involved --------------------
summarise_realdata_rt_estimgamma <- function(stan_posterior,
                                             weekly_data,
                                             start_date,
                                             include_chains){
  
  
  time_week <- weekly_data %>%
    dplyr::select( date, time)
  
  rt_posterior <- stan_posterior %>%
    spread_draws(log_rt[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(time = i + 0 + start_date -1) %>%
    mutate(rt = exp(log_rt)) %>%
    dplyr::select(time, rt) %>%
    group_by(time) %>% 
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    # dplyr::select(date, rt_median = rt, rt_CI95l = .lower, rt_CI95u = .upper) %>%
    left_join(time_week, by = "time")
  
  return(rt_posterior)
}


# create incid posterior for real data ------------------------------------
summarise_realdata_incid_estimgamma <- function(stan_posterior,
                                                weekly_data,
                                                start_date,
                                                include_chains){
  time_week <- weekly_data %>%
    dplyr::select( min_date, time)
  
  incid_posterior <-stan_posterior %>%
    spread_draws(incid[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(i) %>%
    median_qi() %>%
    mutate(date = i + start_date - 1) %>%
    dplyr::select(date, incid_median = incid, incid_CI95l = .lower, incid_CI95u = .upper) %>%
    left_join(time_week, by = c("date" = "time"))
  
  return(incid_posterior)
  
}


# create case posterior predictive for real data --------------------------

summarise_realdata_case_estimgamma <- function(stan_posterior,
                                               weekly_data,
                                               start_date,
                                               include_chains){
  
  true_cases <- weekly_data %>%
    dplyr::select(time, min_date, total_cases)
  
  case_posterior <- stan_posterior %>%
    spread_draws(gen_obs[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(i) %>%
    median_qi() %>%
    mutate(date = i - 0 + start_date - 1) %>%
    dplyr::select(date, case_median = gen_obs, case_CI95l = .lower, case_CI95u = .upper) %>%
    left_join(true_cases, by = c("date" = "time"))
  
  
  return(case_posterior)
  
}

# function for calculating stan rhats by hand -----------------------------

calc_stan_diags <- function(draws, include_chains) {
  draws <- draws %>%
    filter(.chain %in% include_chains) %>%
    dplyr::select(-.draw, 
                  -accept_stat__,
                  -stepsize__,
                  -treedepth__,
                  -n_leapfrog__,
                  -divergent__,
                  -energy__)
  
  num_vars <- dim(draws)[2] - 2
  var_names <- names(draws)[3:dim(draws)[2]]
  rhat <- rep(0, num_vars)
  essbulk <- rep(0, num_vars)
  esstail <- rep(0, num_vars)
  i <- 1
  
  base_columns <- draws %>%
    dplyr::select(.chain, .iteration)
  for (i in 1:num_vars) {
    var_column <- draws[,i+2]
    
    names(var_column) <-  "var"
    
    var_matrix <- as.matrix(cbind(base_columns, var_column) %>%
                              pivot_wider(id_cols = .iteration, 
                                          names_from = .chain,
                                          values_from = var) %>%
                              dplyr::select(-.iteration))
    
    rhat[i] <- as.numeric(Rhat(var_matrix))
    essbulk[i] <- as.numeric(ess_bulk(var_matrix))
    esstail[i] <- as.numeric(ess_tail(var_matrix))
  }
  
  df <- cbind(var_names, rhat, essbulk, esstail)
  range_rhat <- range(rhat)
  range_essbulk <- range(essbulk)
  range_esstail <- range(esstail)
  
  output <- list(df, range_rhat, range_essbulk, range_esstail)
  
  return(output)
  
}

# Discretize  Distributions ------------------------------------------
# Epidemia style discretization of gamma
epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, y)
  pmf[1] <- pgamma(1.5, alpha, rate = beta)
  for (i in 2:y) {
    pmf[i] <- pgamma(i+.5, alpha, rate = beta) - pgamma(i-.5, alpha, rate = beta)
  }
  
  pmf
}

zero_epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, (y+1))
  pmf[1] <- pgamma(0.5, alpha, rate = beta)
  for (i in 2:(y+1)) {
    pmf[i] <- pgamma(i -1 +.5, alpha, rate = beta) - pgamma(i -1 -.5, alpha, rate = beta)
  }
  
  pmf
}


epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, y)
  pmf[1] <- phypoexp(1.5, rates)
  for (i in 2:y) {
    pmf[i] <- phypoexp(i+.5, rates) - phypoexp(i-.5, rates)
  }
  
  pmf
}


epidemia_lognormal <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- plnorm(1.5, meanlog = params[1], sdlog = params[2])
  for (i in 2:y) {
    pmf[i] <- plnorm(i+.5, meanlog = params[1], sdlog = params[2]) - 
      plnorm(i-.5, meanlog = params[1], sdlog = params[2])
  }
  
  pmf
}

epidemia_weibull <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- pweibull(1.5, shape = params[1], scale = params[2])
  for (i in 2:y) {
    pmf[i] <- pweibull(i+.5, shape = params[1], scale = params[2]) - 
      pweibull(i-.5, shape = params[1], scale = params[2])
  }
  
  pmf
}



zero_epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, (y+1))
  pmf[1] <- phypoexp(0.5, rates)
  for (i in 2:y+1) {
    pmf[i] <- phypoexp(i -1 +.5, rates) - phypoexp(i -1 -.5, rates)
  }
  
  pmf
}


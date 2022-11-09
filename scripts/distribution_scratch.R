# distribution scratch


# rho ---------------------------------------------------------------------
log_rho <- rnorm(10000, mean = log(0.2), sd = 0.3)
rho <- exp(log_rho)
quantile(rho, c(0.025, 0.5, 0.975))

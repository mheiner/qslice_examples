# Calculate the log likelihood
dlm_likelihood_components <- function(theta, parameters) {
  # Extract relevant parameters
  ystar <- parameters$ystar
  f <- parameters$f
  m0 <- parameters$m0
  v0 <- parameters$v0
  tau2 <- parameters$tau2
  sigma2 <- parameters$sigma2
  time_diff <- parameters$time_diff
  
  n <- length(ystar)
  
  stopifnot(length(time_diff) == n) # time_diff includes 0 to 1
  stopifnot(length(theta) == n)
  
  stdev_obs <- sqrt(sigma2)

  process_means <- c(m0, theta[1:(n-1)])
  process_sds <- sqrt(c(
    v0 + tau2*time_diff[1], 
    rep(tau2, n-1) * time_diff[2:n]
    )) # assumes theta0 integrated out, time diff from 0 to 1 is 1.0
  lprob_theta_positive <- pnorm(0, process_means, process_sds, lower.tail = FALSE, log.p = TRUE)
  
  llik_proc <- dnorm(theta, mean = process_means, sd = process_sds, log = TRUE) - lprob_theta_positive
  llik_obs <- dnorm(ystar, mean = theta * f, sd = stdev_obs, log = TRUE)
  
  list(llik_proc = llik_proc, llik_obs = llik_obs)
}

dlm_likelihood <- function(theta, parameters) {
  tmp <- dlm_likelihood_components(theta = theta, parameters = parameters)
  sum(tmp$llik_proc) + sum(tmp$llik_obs)
}

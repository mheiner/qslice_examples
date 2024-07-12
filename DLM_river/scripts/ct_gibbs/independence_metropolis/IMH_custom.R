independence_metropolis_update <- function(ystar, f, m0, v0, tau2, sigma2, old_theta, 
                                           time_diff, pseu_family, pseu_params) {
  
  # Put parameters into a list to be supplied to helper functions
  parameters <- list(
    ystar = ystar,
    f = f,
    m0 = m0,
    v0 = v0,
    tau2 = tau2,
    sigma2 = sigma2,
    old_theta = old_theta,
    time_diff = time_diff,
    pseu_family = pseu_family,
    pseu_params = pseu_params
  )
  
  # Extract length of the time series
  n <- length(ystar)
  
  # Calculate Kalman filter
  kalman_filter <- ff(parameters = parameters)
  
  # Draw new proposal
  theta_proposal <- bs_icdf_trunc(runif(n), kalman_filter, pseu_family, pseu_params)

  # Calculate log-likelihood for previous iteration
  old_dlm_log_likelihood <- dlm_likelihood(old_theta, parameters)
  # Calculate the FFBS log-likelihood for previous iteration
  old_ffbs_log_likelihood <- bs_likelihood_constrained(old_theta, kalman_filter, pseu_family, pseu_params)
  # Calculate target log-likelihood
  old_target_log_likelihood <- old_dlm_log_likelihood - old_ffbs_log_likelihood
  
  # Calculate log-likelihood for proposal
  proposal_dlm_log_likelihood <- dlm_likelihood(theta_proposal, parameters)
  # Calculate the FFBS log-likelihood for proposal
  proposal_ffbs_log_likelihood <- bs_likelihood_constrained(theta_proposal, kalman_filter, pseu_family, pseu_params)
  # Calculate target log-likelihood for proposal
  proposal_target_log_likelihood <- proposal_dlm_log_likelihood - proposal_ffbs_log_likelihood
  
  # Calculate log acceptance rate
  log_acceptance <- proposal_target_log_likelihood - old_target_log_likelihood
  
  # Log cutoff draw
  z <- log(runif(1))
  
  # Accept or reject
  accept <- log_acceptance > z
  
  # Return the accepted value or the old value
  if(accept) {
    return(theta_proposal)
  } else {
    return(old_theta)
  }
}


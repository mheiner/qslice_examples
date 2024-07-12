#################################################################
### Slice sampler with backward sampler proposal distribution ###
#################################################################
# "PITS sampler" (probability integral transform slice sampler)

pits_update_constrained <- function(ystar, f, m0, v0, tau2, sigma2, old_theta, time_diff, pseu_family, pseu_params) {
  
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
  
  # Set number of rejections to 0
  rejections <- 0
  
  # Calculate Kalman filter
  kalman_filter <- ff(parameters = parameters)
  # Calculate log-likelihood for previous iteration
  old_dlm_log_likelihood <- dlm_likelihood(old_theta, parameters)
  # Calculate the FFBS log-likelihood for previous iteration
  old_ffbs_log_likelihood <- bs_likelihood_constrained(old_theta, kalman_filter, pseu_family, pseu_params)
  # Calculate target log-likelihood
  old_target_log_likelihood <- old_dlm_log_likelihood - old_ffbs_log_likelihood
  # Sample the standard uniform(0,1)
  z <- runif(1)
  # calculate the log cutoff by log(uniform(0,L)) = log(uniform(0,1) * L) = log(z) + log(L)
  log_cutoff <- log(z) + old_target_log_likelihood
  
  # Create original hyper cube to later be restricted to a hyper rectangle
  hyper_rectangle <- cbind(rep(0, n), rep(1, n))
  # Get old phi vector by applying the CDF to each theta
  old_phi <- bs_cdf_trunc(old_theta, kalman_filter, pseu_family, pseu_params)
  # Sample a proposed phi vector (CDF of theta)
  phi_proposal <- runif(n)
  # Convert proposed phi back to theta
  theta_proposal <- bs_icdf_trunc(phi_proposal, kalman_filter, pseu_family, pseu_params)
  # Check the target likelihood of the proposal
  proposal_dlm_log_likelihood <- dlm_likelihood(theta_proposal, parameters)
  proposal_ffbs_log_likelihood <- bs_likelihood_constrained(theta_proposal, kalman_filter, pseu_family, pseu_params)
  proposal_target_log_likelihood <- proposal_dlm_log_likelihood - proposal_ffbs_log_likelihood
  # Resample the proposal as long as the proposal has lower likelihood than the cutoff
  while((proposal_target_log_likelihood < log_cutoff)) {
    # Increase the number of rejections by 1
    rejections <- rejections + 1
    
    # Restrict the hyper rectangle prior to resampling
    proposal_lower <- phi_proposal < old_phi
    proposal_higher <- !proposal_lower
    if (any(proposal_lower)) {
      hyper_rectangle[proposal_lower,1] <- phi_proposal[proposal_lower]
    }
    if (any(proposal_higher)) {
      hyper_rectangle[proposal_higher,2] <- phi_proposal[proposal_higher]
    }
    # Draw new proposal
    phi_proposal <- runif(n = n, min = hyper_rectangle[,1], max = hyper_rectangle[,2])
    # Convert proposed phi back to theta
    theta_proposal <- bs_icdf_trunc(phi_proposal, kalman_filter, pseu_family, pseu_params)
    # Check the target likelihood of the proposal
    proposal_dlm_log_likelihood <- dlm_likelihood(theta_proposal, parameters)
    proposal_ffbs_log_likelihood <- bs_likelihood_constrained(theta_proposal, kalman_filter, pseu_family, pseu_params)
    proposal_target_log_likelihood <- proposal_dlm_log_likelihood - proposal_ffbs_log_likelihood
  }
  
  # Once a proposal is accepted return it as the final value
  return(
    list(
      theta = theta_proposal,
      phi = phi_proposal,
      rejections = rejections
    )
  )
}

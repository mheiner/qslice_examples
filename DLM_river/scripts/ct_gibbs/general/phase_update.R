phase_update <- function(y, old_phases, sigma2, theta_matrix, frequencies, period) {
  
    # Put parameters into a list to be supplied to helper functions
    parameters <- list(
      y = y,
      sigma2 = sigma2,
      theta_matrix = theta_matrix,
      frequencies = frequencies,
      period = period
    )
  
  # Extract length of the time series
  n <- length(y)
  n_freq <- length(frequencies)

  # Calculate log-likelihood for previous iteration
  old_dlm_log_likelihood <- phase_dlm_likelihood(old_phases, parameters, period)
  # Calculate the FFBS log-likelihood for previous iteration
  old_beta_log_likelihood <- beta_log_likelihood(old_phases, old_phases, period)
  # Calculate target log-likelihood
  old_target_log_likelihood <- old_dlm_log_likelihood - old_beta_log_likelihood
  # Sample the standard uniform(0,1)
  z <- runif(1)
  # calculate the log cutoff by log(uniform(0,L)) = log(uniform(0,1) * L) = log(z) + log(L)
  log_cutoff <- log(z) + old_target_log_likelihood
  
  # Create original hyper cube to later be restricted to a hyper rectangle
  hyper_rectangle <- matrix(c(0,1),
                            nrow = n_freq,
                            ncol = 2,
                            byrow = TRUE)
  # Get old phi vector by applying the CDF to each phase
  old_phi <- beta_cdf(old_phases, old_phases, period)
  # Sample a proposed phi vector (CDF of phase)
  phi_proposal <- runif(n_freq)
  # Convert proposed phi back to phase
  phase_proposal <- beta_quantiles(phi_proposal, old_phases, period)
  # Check the target likelihood of the proposal
  proposal_dlm_log_likelihood <- phase_dlm_likelihood(phase_proposal, parameters, period)
  proposal_beta_log_likelihood <- beta_log_likelihood(phase_proposal, old_phases, period)
  proposal_target_log_likelihood <- proposal_dlm_log_likelihood - proposal_beta_log_likelihood
  # Resample the proposal as long as the proposal has lower likelihood than the cutoff
  while(proposal_target_log_likelihood < log_cutoff) {
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
    phi_proposal <- runif(n = n_freq, min = hyper_rectangle[,1], max = hyper_rectangle[,2])
    # Convert proposed phi back to phase
    phase_proposal <- beta_quantiles(phi_proposal, old_phases, period)
    # Check the target likelihood of the proposal
    proposal_dlm_log_likelihood <- phase_dlm_likelihood(phase_proposal, parameters, period)
    proposal_beta_log_likelihood <- beta_log_likelihood(phase_proposal, old_phases, period)
    proposal_target_log_likelihood <- proposal_dlm_log_likelihood - proposal_beta_log_likelihood
  }
  
  # Once a proposal is accepted return it as the final value
  return(phase_proposal)
}

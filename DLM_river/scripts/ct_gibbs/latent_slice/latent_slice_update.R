#############################
### Slice Sampler for DLM ###
#############################

latent_slice_update <- function(ystar, f, m0, v0, tau2, sigma2, old_s, latent_scale, old_theta, time_diff) {
  
  # Put parameters into a list to be supplied to helper functions
  parameters <- list(
    ystar = ystar,
    f = f,
    m0 = m0,
    v0 = v0,
    tau2 = tau2,
    sigma2 = sigma2,
    old_s = old_s,
    latent_scale = latent_scale,
    old_theta = old_theta,
    time_diff = time_diff
  )
  
  stopifnot(all(old_theta > 0.0))
  
  # Extract length of the time series
  n <- length(ystar)
  
  # Set number of rejections to 0
  rejections <- 0
  
  # Calculate log-likelihood for previous iteration
  old_log_likelihood <- dlm_likelihood(old_theta, parameters)
  # Sample the standard uniform(0,1)
  z <- runif(1)
  # calculate the log cutoff by log(uniform(0,L)) = log(uniform(0,1) * L) = log(z) + log(L)
  log_cutoff <- log(z) + old_log_likelihood
  
  # Halve the s's
  half_s <- old_s/2
  l <- runif(n, old_theta - half_s, old_theta + half_s)
  s <- -log(runif(n)) * latent_scale + 2 * abs(l - old_theta)
  half_s <- s/2
  L <- l - half_s
  U <- l + half_s
  
  L <- pmax(0.0, L) # preempt rejections
  stopifnot(all(L < U))
  
  theta_proposal <- runif(n, L, U)
  proposal_log_likelihood <- dlm_likelihood(theta_proposal, parameters)

  # Resample the proposal as long as the proposal has lower likelihood than the cutoff
  while(proposal_log_likelihood < log_cutoff) {
    
    rejections <- rejections + 1
    
    # Restrict the hyper rectangle prior to resampling
    lower <- theta_proposal < old_theta
    if (any(lower)) {
      L[lower] <- theta_proposal[lower]
    }
    if (any(!lower)) {
      U[!lower] <- theta_proposal[!lower]
    }
    
    theta_proposal <- runif(n, L, U)
    
    proposal_log_likelihood <- dlm_likelihood(theta_proposal, parameters)
    
  }
  
  # Once a proposal is accepted return it as the final value
  return(
    list(
      theta = theta_proposal,
      rejections = rejections,
      s = s
    )
  )
}

#############################
### Slice Sampler for DLM ###
#############################

slice_update <- function(ystar, f, m0, v0, tau2, sigma2, width, old_theta, time_diff) {
  
  # Put parameters into a list to be supplied to helper functions
  parameters <- list(
    ystar = ystar,
    f = f,
    m0 = m0,
    v0 = v0,
    tau2 = tau2,
    sigma2 = sigma2,
    width = width,
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

  ## one option is to set w (width) based on FFBS; didn't work well
  # kfilt <- ff(parameters = parameters)
  # bs_sc <- bs_stdev(theta = old_theta, kalman_filter = kfilt)
  bs_sc <- 1
  
  # Generate a random hyper rectangle around the old theta
  hyper_rectangle <- matrix(nrow = n, ncol = 2)
  # Establish lower bound
  hyper_rectangle[,1] <- old_theta - width * bs_sc * runif(n)
  # Establish upper bound
  hyper_rectangle[,2] <- hyper_rectangle[,1] + width
  # Move lower bound up to zero since we'll reject anyway
  hyper_rectangle[,1] <- pmax(0.0, hyper_rectangle[,1])
  
  # Prepare a numeric vector for the proposed theta
  theta_proposal <- numeric(n)
  # Draw a proposal theta from the hyper rectangle
  theta_proposal <- runif(n = n, min = hyper_rectangle[,1], max = hyper_rectangle[,2])

  # Check the likelihood of the proposal
  proposal_likelihood <- dlm_likelihood(theta_proposal, parameters)
  
  # Resample the proposal as long as the proposal has lower likelihood than the cutoff
  while(proposal_likelihood < log_cutoff) {
    
    rejections <- rejections + 1
    
    # Restrict the hyper rectangle prior to resampling
    proposal_lower <- theta_proposal < old_theta
    proposal_higher <- !proposal_lower
    if (any(proposal_lower)) {
      hyper_rectangle[proposal_lower,1] <- theta_proposal[proposal_lower]
    }
    if (any(proposal_higher)) {
      hyper_rectangle[proposal_higher,2] <- theta_proposal[proposal_higher]
    }
    
    theta_proposal <- runif(n = n, min = hyper_rectangle[,1], max = hyper_rectangle[,2])
    
    proposal_likelihood <- dlm_likelihood(theta_proposal, parameters)
  }
  
  # Once a proposal is accepted return it as the final value
  return(
    list(
      theta = theta_proposal,
      rejections = rejections
    )
  )
}

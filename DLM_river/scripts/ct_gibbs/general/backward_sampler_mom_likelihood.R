bs_moments <- function(theta, kalman_filter) {
  
  ## for one-dimensional theta (at each time point)
  
  # Store filter objects
  filtered_means <- kalman_filter$filtered_means
  filtered_variance <- kalman_filter$filtered_variance
  prior_variance <- kalman_filter$prior_variance
  prior_means <- kalman_filter$prior_mean
  n <- length(filtered_means)
  
  # initialize vectors to hold back sampling information
  # backward sampled mean vector and covariance estimates
  bs_means <- numeric(n)
  bs_variance <- numeric(n)
  
  # Insert the kalman filter values at the end of the vectors
  bs_means[n] <- filtered_means[n]
  bs_variance[n] <- filtered_variance[n]
  
  ### Calculate backward sampler moments: Eqns 4.12, 4.13 in Prado & West (2010)
  bs_variance[1:(n-1)] <- filtered_variance[1:(n-1)] * (1.0 - filtered_variance[1:(n-1)] / prior_variance[2:n])
  bs_means[1:(n-1)] <- filtered_means[1:(n-1)] + (filtered_variance[1:(n-1)] * (theta[2:n] - prior_means[2:n]) / prior_variance[2:n]) # previously had (theta[2:n] - filtered_means[1:(n-1)])
  bs_stdev <- sqrt(bs_variance)
  
  list(means = bs_means, var = bs_variance, stdev = bs_stdev)
}

bs_stdev <- function(theta, kalman_filter) {
  
  # Store filter objects
  filtered_variance <- kalman_filter$filtered_variance
  prior_variance <- kalman_filter$prior_variance
  n <- length(kalman_filter$filtered_means)
  
  # initialize vectors to hold back sampling information
  # backward sampled mean vector and covariance estimates
  bs_variance <- numeric(n)
  
  # Insert the kalman filter values at the end of the vectors
  bs_variance[n] <- filtered_variance[n]
  
  ### Calculate backward sampler moments: Eqns 4.12, 4.13 in Prado & West (2010)
  bs_variance[1:(n-1)] <- filtered_variance[1:(n-1)] * (1.0 - filtered_variance[1:(n-1)] / prior_variance[2:n])

  sqrt(bs_variance)
}

bs_likelihood <- function(theta, kalman_filter) {
  tmp <- bs_moments(theta, kalman_filter)
  
  llik_vec <- dnorm(theta, tmp$means, tmp$stdev, log = TRUE)
  sum(llik_vec)
}

bs_likelihood_constrained <- function(theta, kalman_filter, pseu_family, pseu_params) {
  
  tmp <- bs_moments(theta, kalman_filter)
  
  if (pseu_family == "normal") {
    
    log_normaliz_cnst <- pnorm(0.0, mean = tmp$means, sd = tmp$stdev,
                               lower.tail = FALSE, log.p = TRUE)
    
    llik0 <- dnorm(theta, tmp$means, tmp$stdev, log = TRUE)
    
  } else if (pseu_family == "t") {
    
    require("crch")
    
    log_normaliz_cnst <- ptt(0.0, location = tmp$means, scale = tmp$stdev,
                             df = pseu_params$degf,
                             lower.tail = FALSE, log.p = TRUE)
    
    llik0 <- dtt(theta, location = tmp$means, scale = tmp$stdev, 
                 df = pseu_params$degf,
                 log = TRUE)
    
  }
  
  sum(llik0) - sum(log_normaliz_cnst)
}
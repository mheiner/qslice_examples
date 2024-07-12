bs_qnorm <- function(probabilities, kalman_filter) {
  # Store filter objects
  filtered_means <- kalman_filter$filtered_means
  filtered_variance <- kalman_filter$filtered_variance
  prior_variance <- kalman_filter$prior_variance
  prior_means <- kalman_filter$prior_mean
  n <- length(filtered_means)
  
  # Initialize theta vector
  theta <- numeric(n)
  
  # Get last theta value from final filter distribution
  theta[n] <- qnorm(p = probabilities[n],
                         mean = filtered_means[n], 
                         sd = sqrt(filtered_variance[n]))
  
  # initialize vectors to hold back sampling information
  # backward sampled mean vector and covariance matrix estimates
  bs_means <- numeric(n)
  bs_variance <- numeric(n)
  
  # Insert the kalman filter values at the end of the vectors
  bs_means[n] <- filtered_means[n]
  bs_variance[n] <- filtered_variance[n]
  
  # Calculate the backward sampler variances
  bs_variance[1:(n-1)] <- filtered_variance[1:(n-1)] * (1.0 - filtered_variance[1:(n-1)] / prior_variance[2:n])
  bs_stdev <- sqrt(bs_variance)
  
  ### Calculate backward sampler distributions ###
  
  for(i in (n-1):1) {
    # Using theta from the time point just ahead calculate the conditional mean vector
    bs_means[i] <- filtered_means[i] + (filtered_variance[i] * (theta[i+1] - prior_means[i+1]) / prior_variance[i+1])
    
    # Using the calculated distribution and given probability find theta_i
    theta[i] <- qnorm(p = probabilities[i], a = 0, 
                      mean = bs_means[i], sd = bs_stdev[i])
  }
  
  return(theta)
}

bs_icdf_trunc <- function(probabilities, kalman_filter, pseu_family, pseu_params) {
  # Store filter objects
  filtered_means <- kalman_filter$filtered_means
  filtered_variance <- kalman_filter$filtered_variance
  prior_variance <- kalman_filter$prior_variance
  prior_means <- kalman_filter$prior_mean
  n <- length(filtered_means)
  
  # Initialize theta vector
  theta <- numeric(n)
  
  # Get last theta value from final filter distribution
  if (pseu_family == "normal") {
    theta[n] <- qtruncnorm(p = probabilities[n], a = 0.0, 
                           mean = filtered_means[n], 
                           sd = sqrt(filtered_variance[n]))
  } else if (pseu_family == "t") {
    theta[n] <- crch::qtt(p = probabilities[n], left = 0.0, 
                          location = filtered_means[n], 
                          scale = sqrt(filtered_variance[n]),
                          df = pseu_params$degf)
  }
  
  # initialize vectors to hold back sampling information
  # backward sampled mean vector and covariance matrix estimates
  bs_means <- numeric(n)
  bs_variance <- numeric(n)
  
  # Insert the kalman filter values at the end of the vectors
  bs_means[n] <- filtered_means[n]
  bs_variance[n] <- filtered_variance[n]
  
  # Calculate the backward sampler variances
  bs_variance[1:(n-1)] <- filtered_variance[1:(n-1)] * (1.0 - filtered_variance[1:(n-1)] / prior_variance[2:n])
  bs_stdev <- sqrt(bs_variance)
  
  ### Calculate backward sampler distributions ###
  
  for(i in (n-1):1) {
    # Using theta from the time point just ahead calculate the conditional mean vector
    bs_means[i] <- filtered_means[i] + (filtered_variance[i] * (theta[i+1] - prior_means[i+1]) / prior_variance[i+1])
    
    # Using the calculated distribution and given probability find theta_i
    if (pseu_family == "normal") {
      theta[i] <- qtruncnorm(p = probabilities[i], a = 0, 
                             mean = bs_means[i], sd = bs_stdev[i])      
    } else if (pseu_family == "t") {
      theta[i] <- crch::qtt(p = probabilities[i], left = 0.0, 
                            location = bs_means[i], 
                            scale = bs_stdev[i],
                            df = pseu_params$degf)
    }
  }
  
  return(theta)
}

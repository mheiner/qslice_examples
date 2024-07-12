backward_sample <- function(kalman_filter) {
  # Store filter objects
  filtered_means <- kalman_filter$filtered_means
  filtered_variance <- kalman_filter$filtered_variance
  prior_variance <- kalman_filter$prior_variance
  prior_means <- kalman_filter$prior_mean
  n <- length(filtered_means)
  
  ###
  ### Back sample to get estimates based on all data
  ###
  
  ### initialize vectors to hold back sampled information ###
  # backs ampled mean vector and covariance matrix estimates
  mean_back_sample <- numeric(n)
  var_back_sample <- numeric(n)
  
  ### Insert the kalman filter values at the end of the vector ###
  mean_back_sample[n] <- filtered_means[n]
  var_back_sample[n] <- filtered_variance[n]
  
  var_back_sample[1:(n-1)] <- filtered_variance[1:(n-1)] * (1 - filtered_variance[1:(n-1)] / prior_variance[2:n])
  sd_back_sample <- sqrt(var_back_sample)
  
  ### Back sample once ###
  theta <- numeric(n)
  # Generate a possible theta vector from the final filtered distribution
  theta[n] <- sd_back_sample[n] * rnorm(1) + filtered_means[n]
  
  for(i in (n-1):1) {
    # Using the draw from the time point just ahead calculate the conditional mean vector
    mean_back_sample[i] <- filtered_means[i] + filtered_variance[i] * (theta[i + 1] - prior_means[i + 1]) / prior_variance[i + 1]
    # Draw from the distribution created by the conditional mean and smoothed variance
    theta[i] <- sd_back_sample[i] * rnorm(1) + mean_back_sample[i]
  }
  
  return(theta)
}
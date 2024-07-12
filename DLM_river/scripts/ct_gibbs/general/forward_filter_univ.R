ff <- function(ystar, f, m0, v0, tau2, sigma2, time_diff, parameters = NULL) {
  
  if (!is.null(parameters)) {
    # Extract relevant parameters
    ystar <- parameters$ystar
    f <- parameters$f
    m0 <- parameters$m0
    v0 <- parameters$v0
    tau2 <- parameters$tau2
    sigma2 <- parameters$sigma2
    time_diff <- parameters$time_diff
  }
  
  ###
  ### Calculate distribution of yt|y1-yt for all t (Kalman Filter)
  ###
  
  # Calculate useful quantities from input
  n <- length(ystar)
  stopifnot(length(time_diff) == n)
  
  ### Initialize vectors to hold Kalman filter information ###
  # Prior mean vectors and covariance matrices
  mean_prior <- numeric(n)
  var_prior <- numeric(n)
  # Posterior mean vectors and covariance matrices
  mean_filter <- numeric(n)
  var_filter <- numeric(n)
  # Predictive distribution values and adaptive vector
  mean_pred <- numeric(n)
  var_pred <- numeric(n)
  adaptive <- numeric(n)
  
  ### Start by calculating the first distribution ###
  # Prior for the theta vector
  mean_prior[1] <- m0
  var_prior[1] <- v0 + tau2 * time_diff[1]
  # Prediction mean and variance for the next observation
  mean_pred[1] <- f[1] * mean_prior[1]
  var_pred[1] <- (f[1])^2 * var_prior[1] + sigma2
  # Adaptive matrix to simplify filter calculations (Prado/West 125)
  adaptive[1] <- (var_prior[1] * f[1]) / var_pred[1]
  # Mean vector and covariance matrix for the filter
  mean_filter[1] <- mean_prior[1] + adaptive[1] * (ystar[1] - mean_pred[1])
  var_filter[1] <- var_prior[1] - (adaptive[1])^2 * var_pred[1]
  
  ### Now calculate for all other t ###
  for(i in 2:n) {
    # Generate the prior for the next mean vector and covariance matrix
    mean_prior[i] <- mean_filter[i - 1]
    var_prior[i] <- var_filter[i - 1] + tau2 * time_diff[i]
    # Mean and variance for the next observation
    mean_pred[i] <- f[i] * mean_prior[i]
    var_pred[i] <- (f[i])^2 * var_prior[i] + sigma2
    # Adaptive matrix for the next observation
    adaptive[i] <- (var_prior[i] * f[i]) / var_pred[i]
    # Filtered mean vector and covariance matrix
    mean_filter[i] <- mean_prior[i] + adaptive[i] * (ystar[i] - mean_pred[i])
    var_filter[i] <- var_prior[i] - (adaptive[i])^2 * var_pred[i]
  }
  
  return(
    list(
      filtered_means = mean_filter,
      filtered_variance = var_filter,
      prior_variance = var_prior,
      prior_mean = mean_prior
    )
  )
}

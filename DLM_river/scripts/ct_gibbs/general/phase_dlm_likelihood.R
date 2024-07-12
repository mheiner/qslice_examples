# Calculate the log likelihood
phase_dlm_likelihood <- function(phases, parameters, period) {
  # Extract relevant parameters
  y <- parameters$y
  sigma2 <- parameters$sigma2
  theta_matrix <- parameters$theta_matrix
  frequencies <- parameters$frequencies
  period <- parameters$period
  
  # Convert the phases to F_mat
  n_freq <- length(phases)
  n <- length(y)
  F_mat <- matrix(nrow = n, ncol = n_freq + 1)
  F_mat[,1] <- 1
  F_mat[,2:(n_freq + 1)] <- cos(2 * pi * t(t(times %o% frequencies) + phases) / period)

  # Calculate the observation means by multiplying each theta at each time point by
  #  its corresponding F and summing
  observation_means <- colSums(t(F_mat) * theta_matrix)
  # Calculate the log likelihood
  log_likelihood <- sum(dnorm(y, mean = observation_means, sd = sqrt(sigma2), log = TRUE))
  
  # Return the log-likelihood
  return(log_likelihood)
}

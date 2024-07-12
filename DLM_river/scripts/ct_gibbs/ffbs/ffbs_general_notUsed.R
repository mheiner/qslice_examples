ffbs <- function(y, F_mat, G, m0, C0, W, sigma2, time_diff) {
  ###
  ### Calculate E(yt|y1-yt) for all t (Kalman Filter)
  ###
  
  # Calculate useful quantities from input
  n <- length(y)
  p <- ncol(F_mat)
  
  # Scale W by the time differences
  Wt <- array(W, dim = c(dim(W),n))
  for(i in 1:n) Wt[,,i] <- Wt[,,i] * time_diff[i]
  
  ### Initialize vectors to hold Kalman filter information ###
  # Prior mean vectors and covariance matrices
  mean_prior <- matrix(NA, nrow = n, ncol = p)
  var_prior <- vector("list",n)
  var_prior_chol <- vector("list",n)
  # Posterior mean vectors and covariance matrices
  mean_filter <- matrix(NA, nrow = n, ncol = p)
  var_filter <- vector("list",n)
  var_filter_chol <- vector("list",n)
  # Predictive distribution values and adaptive vector
  mean_pred <- numeric(n)
  var_pred <- numeric(n)
  adaptive_mat <- matrix(NA, nrow = n, ncol = p)
  
  ### Start by calculating the first distribution ###
  # Prior for the next theta vector
  mean_prior[1,] <- G %*% m0
  var_prior[[1]] <- G %*% C0 %*% t(G) + Wt[,,1]
  var_prior_chol[[1]] <- chol(var_prior[[1]])
  # Prediction mean and variance for the next observation
  mean_pred[1] <- crossprod(F_mat[1,], mean_prior[1,])
  var_pred[1] <- tcrossprod(tcrossprod(F_mat[1,],var_prior_chol[[1]])) + sigma2
  # Adaptive matrix to simplify filter calculations (Prado/West 125)
  adaptive_mat[1,] <- (var_prior[[1]] %*% F_mat[1,])/var_pred[1]
  # Mean vector and covariance matrix for the filter
  mean_filter[1,] <- mean_prior[1,] + adaptive_mat[1,] * (y[1] - mean_pred[1])
  var_filter[[1]] <- var_prior[[1]] - tcrossprod(adaptive_mat[1,]) * var_pred[1]
  var_filter_chol[[1]] <- chol(var_filter[[1]])
  
  ### Now calculate for all other t ###
  for(i in 2:n) {
    # Generate the prior for the next mean vector and covariance matrix
    mean_prior[i,] <- G %*% mean_filter[i-1,]
    var_prior[[i]] <- tcrossprod(tcrossprod(G,var_filter_chol[[i-1]])) + Wt[,,i]
    var_prior_chol[[i]] <- chol(var_prior[[i]])
    # Mean and variance for the next observation
    mean_pred[i] <- crossprod(F_mat[i,], mean_prior[i,])
    var_pred[i] <- tcrossprod(tcrossprod(F_mat[i,], var_prior_chol[[i]])) + sigma2
    # Adaptive matrix for the next observation
    adaptive_mat[i,] <- (var_prior[[i]] %*% F_mat[i,])/var_pred[i]
    # Filtered mean vector and covariance matrix
    mean_filter[i,] <- mean_prior[i,] + adaptive_mat[i,] * (y[i] - mean_pred[i])
    var_filter[[i]] <- var_prior[[i]] - tcrossprod(adaptive_mat[i,]) * var_pred[i]
    var_filter_chol[[i]] <- chol(var_filter[[i]])
  }
  
  ###
  ### Back sample to get estimates based on all data
  ###
  
  ### initialize vectors to hold smoothing information ###
  # Smoothed mean vector and covariance matrix estimates
  mean_smoothed <- matrix(nrow = n, ncol = p)
  var_smoothed <- vector("list", n)
  var_smoothed_chol <- vector("list", n)
  
  ### Insert the kalman filter values at the end of the matrix/list ###
  mean_smoothed[n,] <- mean_filter[n,]
  var_smoothed[[n]] <- var_filter[[n]]
  var_smoothed_chol[[n]] <- var_filter_chol[[n]]
  
  ### Back sample once ###
  theta_matrix <- matrix(nrow = p, ncol = n)
  # Generate a possible theta vector from the final filtered distribution
  theta_matrix[,n] <- crossprod(var_filter_chol[[n]], rnorm(p)) + mean_filter[n,]
  
  for(i in (n-1):1) {
    # Calculate back-sampling variance at time i
    # FIXME See if forwardsolve and backsolve are faster or more stable (double check how they interact with the upper_tri argument and the original setup)
    var_smoothed[[i]] <- var_filter[[i]] - crossprod(G %*% var_filter[[i]], solve(var_prior[[i + 1]], G %*% var_filter[[i]]))
    var_smoothed_chol[[i]] <- chol(var_smoothed[[i]])
    
    # Using the draw from the time point just ahead calculate the conditional mean vector
    mean_smoothed[i,] <- mean_filter[i,] + crossprod(G %*% var_filter[[i]], solve(var_prior[[i + 1]], theta_matrix[,i+1] - G %*% mean_filter[i,]))
    # Draw from the distribution created by the conditional mean and smoothed variance
    theta_matrix[,i] <- crossprod(var_smoothed_chol[[i]], rnorm(p)) + mean_smoothed[i,]
  }
  return(theta_matrix)
}

DLM_sim <- function(C0,m0,v,F_mat,W,G) {
  # Extract other useful values
  n <- nrow(F_mat)
  p <- ncol(F_mat)
  
  # Initialize state and observation objects
  theta <- matrix(NA, nrow = n, ncol = p)
  y <- numeric(n)
  
  ### Draw initial theta and y ###
  # Multivariate normal draw for the theta vector
  theta[1,] <- crossprod(chol(C0), rnorm(p)) + m0
  # Conversion of the drawn theta vector to an observation using F and V values
  # Define observation sd once
  observation_sd <- sqrt(v)
  y[1] <- crossprod(F_mat[1,], theta[1,]) + rnorm(1,0,observation_sd)
  
  # Draw all subsequent thetas and y's
  # Define inovation cholesky decomposition once
  innovation_cholesky <- chol(W)
  for(i in 2:n) {
    # Generate next theta by passing through evolution matrix and adding innovation variance
    theta[i,] <- G %*% theta[i-1,] + crossprod(innovation_cholesky, rnorm(p))
    # Convert thetas to y's using F and V values
    y[i] <- crossprod(F_mat[i,], theta[i,]) + rnorm(1,0,observation_sd)
  }
  
  ### Plot thetas and y's on the same axis ###
  # Initialize state vector
  state <- numeric(n)
  for(i in 1:n) {
    # Condense theta vector into a state in response space
    state[i] <- crossprod(F_mat[i,], theta[i,])
  }
  return(list(true_mean = state, theta = theta, obs = y))
}

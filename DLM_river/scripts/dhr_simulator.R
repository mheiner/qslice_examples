simulate_dhr <- function(n, frequencies, phases, period, C0, m0, v, W, continuous_time = FALSE, parameters = NULL) {
  
  # Handle potential list input
  if(!is.null(parameters)) {
    env <- environment()
    list2env(parameters, envir = env)
  }
  
  if(continuous_time) {
    times <- 1:n + runif(n, -.5,.5)
    time_diff <- c(1, diff(times))
  }
  
  # Extract other useful values
  n_freq <- length(frequencies)
  p <- n_freq + 1
  
  # Initialize state and observation objects
  theta <- matrix(NA, nrow = p, ncol = n)
  y <- numeric(n)
  state <- numeric(n)
  
  # Create the F matrix
  F_mat <- matrix(nrow = n, ncol = p)
  F_mat[,1] <- 1
  if(continuous_time) {
    for(i in 2:p) {
      F_mat[,i] <- cos(2*pi*(times*frequencies[i-1] + phases[i-1])/period)
    }
  } else {
    for(i in 2:p) {
      F_mat[,i] <- cos(2*pi*((1:n)*frequencies[i-1] + phases[i-1])/period)
    }    
  }
  
  
  ### Draw initial theta and y ###
  # draw for theta_1
  if(continuous_time) {
    theta[1,1] <- m0[1] + rnorm(1, sd = sqrt(diag(C0)[1] * time_diff[1]))
    for(i in 2:p){
      draw <- m0[i] + rnorm(1, sd = sqrt(diag(C0)[i] * time_diff[1]))
      while(draw < 0) {
        draw <- m0[i] + rnorm(1, sd = sqrt(diag(C0)[i] * time_diff[1]))
      }
      theta[i,1] <- draw
    }
  } else {
    theta[1,1] <- m0[1] + rnorm(1, sd = sqrt(diag(C0)[1]))
    for(i in 2:p){
      draw <- m0[i] + rnorm(1, sd = sqrt(diag(C0)[i]))
      while(draw < 0) {
        draw <- m0[i] + rnorm(1, sd = sqrt(diag(C0)[i]))
      }
      theta[i,1] <- draw
    }
  }
  
  
  
  # Conversion of the drawn theta vector to an observation using F and V values
  # Define observation sd once
  observation_sd <- sqrt(v)
  state[1] <- crossprod(F_mat[1,], theta[,1])
  y[1] <- state[1] + rnorm(1,0,observation_sd)
  
  # Draw all subsequent thetas and y's
  if(continuous_time) {
    for(i in 2:n) {
      # Generate next theta by passing through evolution matrix and adding innovation variance
      theta[1, i] <- theta[1, i-1] + rnorm(1, sd = sqrt(diag(W)[1] * time_diff[i]))
      for(j in 2:p){
        draw <- theta[j,i-1] + rnorm(1, sd = sqrt(diag(W)[j] * time_diff[i]))
        while(draw < 0) {
          draw <- theta[j,i-1] + rnorm(1, sd = sqrt(diag(W)[j] * time_diff[i]))
        }
        theta[j,i] <- draw
      }
      # Convert thetas to y's using F and V values
      state[i] <- crossprod(F_mat[i,], theta[,i])
      y[i] <- state[i] + rnorm(1,0,observation_sd)
    }
  } else {
    for(i in 2:n) {
      # Generate next theta by passing through evolution matrix and adding innovation variance
      theta[1,i] <- theta[1,i-1] + rnorm(1, sd = sqrt(diag(W)[1]))
      for(j in 2:p){
        draw <- theta[j,i-1] + rnorm(1, sd = sqrt(diag(W)[j]))
        while(draw < 0) {
          draw <- theta[j,i-1] + rnorm(1, sd = sqrt(diag(W)[j]))
        }
        theta[j,i] <- draw
      }
      # Convert thetas to y's using F and V values
      state[i] <- crossprod(F_mat[i,], theta[,i])
      y[i] <- state[i] + rnorm(1,0,observation_sd)
    }
  }
  
  if(continuous_time) {
    return(
      list(
        true_mean = state,
        theta = theta,
        obs = y,
        times = times
      )
    )
  } else {
    return(
      list(
        true_mean = state,
        theta = theta,
        obs = y
      )
    )
  }
}

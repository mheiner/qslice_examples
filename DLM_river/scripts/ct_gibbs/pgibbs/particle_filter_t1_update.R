###############################
### Particle Filter for DLM ###
###############################

particle_filter_t1_update <- function(n_particles, ystar, f, m0, v0, tau2, sigma2, old_theta, time_diff) {

  require("truncnorm")
    
  # Extract length of the time series
  n <- length(ystar)
  stopifnot(length(time_diff) == n)

  sd_obs <- sqrt(sigma2)
  sd_sys <- sqrt(c(
    v0 + tau2*time_diff[1],
    tau2 * time_diff[2:n] 
  )) # already includes time diff from 0 to 1

  # Create a matrix to store the particles and a matrix to track ancestry
  Np1 <- n_particles + 1
  particles <- matrix(nrow = Np1, ncol = n)
  weights <- matrix(1.0, nrow = Np1, ncol = n)
  A <- matrix(1, nrow = Np1, ncol = n) # allocations for resampling (first column will never be used, but we keep it for easy indexing)
  
  # Generate the particle system
  # Generate the starting particles
  particles[1, 1:n] <- old_theta
  particles[2:Np1, 1] <- rtruncnorm(n_particles, mean = m0, sd = sd_sys[1], a = 0)
  weights[,1] <- dnorm(ystar[1], mean = particles[,1] * f[1], sd = sd_obs)

  ## Bootstrap filter with importance resampling at every step
  ## as implemented in Algo 16.5 of Chopin
  ## the star trajectory is the old theta vector (per Algo 16.7, the wider Gibbs sampler)
  
  for(tt in 2:n) {
    
    A[2:Np1, tt] <- sample.int(Np1, size = n_particles, replace = TRUE, prob = weights[,tt-1]) # fixes a typo(?) in Algo 15.5
    
    # Generate proposed particle trajectory to next time point
    particles[2:Np1, tt] <- rtruncnorm(n_particles, mean = particles[A[2:Np1, tt], tt-1], sd = sd_sys[tt], a = 0) # when tt = 2 this is propagating from time 1 to 2
    
    # Weight the new particles by their likelihoods
    weights[, tt] <- dnorm(ystar[tt], mean = particles[,tt] * f[tt], sd = sd_obs)
  } # End the particle system loop
  
  # Backwards sampling for star trajectory (see Algorithm 16.8 of Chopin)
  star_trajectory <- numeric(n)
  B <- numeric(n)
  B[n] <- sample.int(Np1, size = 1, prob = weights[,n])
  star_trajectory[n] <- particles[B[n], n]
  for(tt in (n-1):1) {
    state_likelihoods <- dtruncnorm(star_trajectory[tt+1], mean = particles[,tt], sd = sd_sys[tt+1], a = 0)
    B[tt] <- sample.int(Np1, size = 1, prob = weights[,tt] * state_likelihoods)
    star_trajectory[tt] <- particles[B[tt], tt]
  }

  # Once a proposal is accepted return it as the final value
  return(star_trajectory)
}

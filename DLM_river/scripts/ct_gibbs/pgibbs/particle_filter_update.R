###############################
### Particle Filter for DLM ###
###############################

particle_filter_update <- function(
  n_particles,
  ystar,
  f,
  m0,
  v0,
  tau2,
  sigma2,
  old_theta,
  theta0,
  time_diff
) {
  # Extract length of the time series
  n <- length(ystar)
  stopifnot(length(time_diff) == n)

  # Create a matrix to store the particles and a matrix to track ancestry
  Np1 <- n_particles + 1
  particles <- matrix(nrow = Np1, ncol = n + 1)
  logweights <- matrix(1.0, nrow = Np1, ncol = n + 1) # weights at time 0 are all proportional to 1
  A <- matrix(1, nrow = Np1, ncol = n) # allocations for resampling

  # Generate the particle system
  # Generate the starting particles
  particles[1, 1] <- theta0 # this is time 0
  particles[2:Np1, 1] <- truncnorm::rtruncnorm(
    n_particles,
    mean = m0,
    sd = sqrt(v0),
    a = 0
  )
  particles[1, 2:(n + 1)] <- old_theta

  ## Bootstrap filter with importance resampling at every step
  ## as implemented in Algo 16.5 of Chopin
  ## the star trajectory is the old theta vector (per Algo 16.7, the wider Gibbs sampler)

  sd_obs <- sqrt(sigma2)
  sd_sys <- sqrt(tau2 * time_diff) # already includes time diff from 0 to 1

  for (tt in 1:n) {
    A[2:Np1, tt] <- sample.int(
      Np1,
      size = n_particles,
      replace = TRUE,
      prob = exp(logweights[, tt] - max(logweights[, tt]))
    ) # fixes a typo(?) in Algo 15.5

    # Generate proposed particle trajectory to next time point
    particles[2:Np1, tt + 1] <- truncnorm::rtruncnorm(
      n_particles,
      mean = particles[A[2:Np1, tt], tt],
      sd = sd_sys[tt],
      a = 0
    ) # when tt = 1 this is propagating from time 0 to time 1

    # Weight the new particles by their likelihoods
    logweights[, tt + 1] <- dnorm(
      ystar[tt],
      mean = particles[, tt + 1] * f[tt],
      sd = sd_obs,
      log = TRUE
    )
  } # End the particle system loop

  # Backwards sampling for star trajectory (see Algorithm 16.8 of Chopin)
  star_trajectory <- numeric(n)
  B <- numeric(n)
  B[n] <- sample.int(
    Np1,
    size = 1,
    prob = exp(logweights[, n + 1] - max(logweights[, n + 1]))
  )
  star_trajectory[n] <- particles[B[n], n + 1]
  for (tt in (n - 1):1) {
    state_log_likelihoods <- dnorm(
      star_trajectory[tt + 1],
      mean = particles[, tt + 1],
      sd = sd_sys[tt + 1],
      log = TRUE
    ) -
      pnorm(
        0,
        mean = particles[, tt + 1],
        sd = sd_sys[tt + 1],
        lower.tail = FALSE,
        log.p = TRUE
      ) # particles has n+1 col, so tt+1 col corresponds to time tt; ALGORITHM AS WRITTEN
    lp_now <- logweights[, tt + 1] + state_log_likelihoods
    B[tt] <- sample.int(
      Np1,
      size = 1,
      prob = exp(lp_now - max(lp_now))
    )
    star_trajectory[tt] <- particles[B[tt], tt + 1]
  }
  state_log_likelihoods <- dnorm(
    star_trajectory[1],
    mean = particles[, 1],
    sd = sd_sys[1],
    log = TRUE
  ) -
    pnorm(
      0,
      mean = particles[, 1],
      sd = sd_sys[1],
      lower.tail = FALSE,
      log.p = TRUE
    )
  b0 <- sample.int(
    Np1,
    size = 1,
    prob = exp(state_log_likelihoods - max(state_log_likelihoods))
  ) # weights[,1] at time 0 is all 1s
  theta0 <- particles[b0, 1]

  # Once a proposal is accepted return it as the final value
  return(
    list(theta = star_trajectory, theta0 = theta0)
  )
}

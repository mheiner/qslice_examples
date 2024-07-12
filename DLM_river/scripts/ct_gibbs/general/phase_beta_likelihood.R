beta_log_likelihood <- function(phases, old_phases, period) {
  
  # Write a function to evaluate probability based on reflection point
  evaluator <- function (phase, old_phase, period) {
    
    # Scale the phase
    scaled_phase <- phase/period
    scaled_old_phase <- old_phase/period
    
    # Find the reflection point,
    # the point at which the edge of the distribution has been shifted to
    if (scaled_old_phase < 0.5) {
      
      scaled_reflection_point <- scaled_old_phase + 0.5
      
    } else {
      
      scaled_reflection_point <- scaled_old_phase - 0.5
      
    }
    
    if (scaled_phase < scaled_reflection_point) {
      
      log_likelihood <- dbeta(scaled_phase + 1 - scaled_reflection_point, 2, 2, log = TRUE)
      
    } else {
      
      log_likelihood <- dbeta(scaled_phase - scaled_reflection_point, 2, 2, log = TRUE)
      
    }
    
    return(log_likelihood)
  }
  
  # Loop through phases and old phases to calculate each likelihood
  n <- length(phases)
  log_likelihoods <- numeric(n)
  for (i in 1:n) {
    log_likelihoods[i] <- evaluator(phases[i], old_phases[i], period)
  }
  
  log_likelihood <- sum(log_likelihoods)
  
  return(log_likelihood)
}

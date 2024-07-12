beta_cdf <- function(phases, old_phases, period) {
  
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
      
      probability <- pbeta(scaled_phase + 1 - scaled_reflection_point, 2, 2)
      
    } else {
      
      probability <- pbeta(scaled_phase - scaled_reflection_point, 2, 2)
      
    }
    
    return(probability)
  }
  
  # Loop through phases and old phases to calculate each probability
  n <- length(phases)
  probabilities <- numeric(n)
  for (i in 1:n) {
    probabilities[i] <- evaluator(phases[i], old_phases[i], period)
  }
  
  return(probabilities)
}

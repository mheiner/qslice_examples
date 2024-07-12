beta_quantiles <- function(probabilities, old_phases, period) {
  
  # Write a function to evaluate probability based on reflection point
  evaluator <- function (probability, old_phase, period) {
    
    # Scale the old phase
    scaled_old_phase <- old_phase/period
    
    # Calculate raw quantile
    raw_quantile <- qbeta(probability, 2, 2)
    
    # Return to original scale
    centered_quantile <- raw_quantile - 0.5 + scaled_old_phase
    if (centered_quantile < 0) {
      constrained_quantile <- centered_quantile + 1
    } else if (centered_quantile > 1) {
      constrained_quantile <- centered_quantile - 1
    } else {
      constrained_quantile <- centered_quantile
    }
    
    quantile <- constrained_quantile * period
    
    return(quantile)
  }
  
  # Loop through probabilities and old phases to calculate each probability
  n <- length(probabilities)
  quantiles <- numeric(n)
  for (i in 1:n) {
    quantiles[i] <- evaluator(probabilities[i], old_phases[i], period)
  }
  
  return(quantiles)
}

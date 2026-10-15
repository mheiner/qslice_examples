bs_pnorm <- function(theta, kalman_filter) {
  tmp <- bs_moments(theta, kalman_filter)

  probabilities <- pnorm(theta, tmp$means, tmp$stdev)

  return(probabilities)
}


bs_pnorm_trunc <- function(theta, kalman_filter) {
  tmp <- bs_moments(theta, kalman_filter)

  probabilities <- .ptnorm0(q = theta, location = tmp$means, scale = tmp$stdev)

  return(probabilities)
}

bs_pt_trunc <- function(theta, kalman_filter, df) {
  tmp <- bs_moments(theta, kalman_filter)

  probabilities <- .ptt0(
    q = theta,
    location = tmp$means,
    scale = tmp$stdev,
    df = df
  )

  return(probabilities)
}

bs_cdf_trunc <- function(theta, kalman_filter, pseu_family, pseu_params) {
  if (pseu_family == "normal") {
    out <- bs_pnorm_trunc(theta = theta, kalman_filter = kalman_filter)
  } else if (pseu_family == "t") {
    out <- bs_pt_trunc(
      theta = theta,
      kalman_filter = kalman_filter,
      df = pseu_params$df
    )
  }

  out
}

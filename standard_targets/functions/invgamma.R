# functions to evalue the inverse gamma distribution

pinvgamma <- function(q, shape, scale) {
  ifelse(q <= 0.0, 0.0, 1 - pgamma(1.0 / q, shape = shape, rate = scale))
}

qinvgamma <- function(p, shape, scale) {
  1.0 / qgamma(1.0 - p, shape = shape, rate = scale)
}

dinvgamma <- function(x, shape, scale, log = FALSE) {
  sapply(x, \(x) {
    if (log) {
      if (x <= 0.0) {
        return(-Inf)
      }
      out <- dgamma(1.0 / x, shape = shape, rate = scale, log = TRUE) -
        2.0 * log(x)
    } else {
      if (x <= 0.0) {
        return(0.0)
      }
      out <- dgamma(1.0 / x, shape = shape, rate = scale, log = FALSE) / x^2
    }
    out
  })
}


rinvgamma <- function(n, shape, scale) {
  1.0 / rgamma(n, shape = shape, rate = scale)
}


psqrtinvgamma <- function(q, igsh, igsc) {
  ifelse(q <= 0.0, 0.0, pinvgamma(q^2, shape = igsh, scale = igsc))
}

qsqrtinvgamma <- function(p, igsh, igsc) {
  qinvgamma(p, shape = igsh, scale = igsc) |> sqrt()
}

dsqrtinvgamma <- function(x, igsh, igsc, log = FALSE) {
  sapply(x, \(x) {
    if (log) {
      if (x <= 0.0) {
        out <- -Inf
      } else {
        out <- dinvgamma(x^2, shape = igsh, scale = igsc, log = TRUE) +
          log(2.0) +
          log(x)
      }
    } else {
      if (x <= 0.0) {
        out <- 0.0
      } else {
        out <- dinvgamma(x^2, shape = igsh, scale = igsc, log = FALSE) * 2.0 * x
      }
    }
  })
  out
}

## authors: Sam Johnnson, Matt Heiner

if (n_extra > 0) {
  ld_use <- function(x) {
    doextra(n_extra)
    dnorm(x, log = TRUE)
  }
} else {
  ld_use <- function(x) {
    dnorm(x, log = TRUE)
  }
}

truth <- list(
  d = function(x) dnorm(x),
  ld = ld_use,
  p = function(x) pnorm(x),
  q = function(u) qnorm(u),
  lb = -Inf,
  ub = Inf,
  t = "normal(0,1)"
)

xlim_range <- c(-4, 4)
ylim_range <- c(0, 0.42)

pseudo_init <- function() c(0.5, 2.0)

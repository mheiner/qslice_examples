## authors: Sam Johnnson, Matt Heiner

truth <- list(
  d = function(x) dnorm(x),
  ld = function(x) dnorm(x, log = TRUE),
  p = function(x) pnorm(x),
  q = function(u) qnorm(u),
  lb = -Inf,
  ub = Inf,
  t = "normal(0,1)"
)

xlim_range <- c(-4, 4)
ylim_range <- c(0, 0.42)

pseudo_init <- function() c(0.5, 2.0)

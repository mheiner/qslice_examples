gshape <- 2.0

truth <- list(
  d = function(x) {
    dgamma(exp(-x), shape = gshape) * exp(-x)
  },
  ld = function(x) {
    dgamma(exp(-x), shape = gshape, log = TRUE) - x
  },
  p = function(x) {
    1.0 - pgamma(exp(-x), shape = gshape)
  },
  q = function(u) {
    -log(qgamma(1.0 - u, shape = gshape))
  },
  lb = -Inf,
  ub = Inf,
  t = paste0("inv-gamma_log(shape = ", gshape, ")")
)

pseudo_init <- function() c(0.5, 2.0)

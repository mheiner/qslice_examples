gshape <- 2.5

if (n_extra > 0) {
  ld_use <- function(x) {
    doextra(n_extra)
    dgamma(exp(x), gshape, log = TRUE) + x
  }
} else {
  ld_use <- function(x) {
    dgamma(exp(x), gshape, log = TRUE) + x
  }
}

truth <- list(
  d = function(x) {
    dgamma(exp(x), gshape) * exp(x)
  },
  ld = ld_use,
  p = function(x) {
    pgamma(exp(x), shape = gshape)
  },
  q = function(u) {
    log(qgamma(u, shape = gshape))
  },
  lb = -Inf,
  ub = Inf,
  t = paste0("gamma_log(shape = ", gshape, ")")
)

pseudo_init <- function() c(0.5, 2.0)

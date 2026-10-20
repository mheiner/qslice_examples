## Set up for the cauch transform
## author: Sam Johnnson

gshape <- 2.5

if (n_extra > 0) {
  ld_use <- function(x) {
    doextra(n_extra)
    dgamma(x, shape = gshape, scale = 1.0, log = TRUE) +
      ifelse(x > 0.0, 0.0, -Inf)
  }
} else {
  ld_use <- function(x) {
    dgamma(x, shape = gshape, scale = 1.0, log = TRUE) +
      ifelse(x > 0.0, 0.0, -Inf)
  }
}

truth <- list(
  d = function(x) {
    dgamma(x, shape = gshape, scale = 1.0) * (x > 0.0)
  },
  ld = ld_use,
  p = function(x) {
    pgamma(x, shape = gshape, scale = 1.0)
  },
  q = function(u) {
    qgamma(u, shape = gshape, scale = 1.0)
  },
  t = paste0("gamma(shape = ", gshape, ")"),
  lb = 0.0,
  ub = Inf
)

xlim_range <- c(-4, 15)
ylim_range <- c(0, 0.42)

pseudo_init <- function() c(0.5, 3.0)

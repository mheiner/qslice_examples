# Tests for stable_decreasing_beta_helpers.R

source("functions/functions_truncbeta.R")
source("functions/functions_QSS.R")

## Sequential U -> V -> U round-trip tests.
set.seed(42)
n_tests <- 1000L
shape1 <- c(2.5, 3.0, 4.0, 2.0)
shape2 <- c(8.0, 6.0, 3.0, 1.5)

round_trip_errors <- replicate(n_tests, {
  u <- runif(length(shape1))
  transformed <- V_decseq_from_U(
    u_vec = u,
    shape1_vec = shape1,
    shape2_vec = shape2
  )
  recovered <- U_from_Vdecseq(
    v_vec = transformed$v,
    shape1_vec = shape1,
    shape2_vec = shape2
  )
  max(abs(recovered$u - u))
})

max_round_trip_error <- max(round_trip_errors)
stopifnot(max_round_trip_error < 1e-10)

## Extreme upper-tail interval where direct lower-tail subtraction is zero.
lower <- 0.999999998
upper <- 0.999999999
shape1_extreme <- 2.0
shape2_extreme <- 50.0

ordinary_interval_probability <-
  pbeta(upper, shape1_extreme, shape2_extreme) -
  pbeta(lower, shape1_extreme, shape2_extreme)

stable_quantile <- qtruncbeta(
  p = 0.5,
  shape1 = shape1_extreme,
  shape2 = shape2_extreme,
  lower = lower,
  upper = upper
)
stable_cdf <- ptruncbeta(
  q = stable_quantile$q,
  shape1 = shape1_extreme,
  shape2 = shape2_extreme,
  lower = lower,
  upper = upper
)

stopifnot(
  ordinary_interval_probability == 0.0,
  is.finite(stable_quantile$log_prob_interval),
  abs(stable_cdf$p - 0.5) < 1e-6
)

cat(
  "Sequential round-trip tests:",
  n_tests,
  "\n",
  "Maximum U -> V -> U error:",
  format(max_round_trip_error, digits = 16),
  "\n",
  "Direct extreme-tail interval probability:",
  ordinary_interval_probability,
  "\n",
  "Stable extreme-tail log interval probability:",
  stable_quantile$log_prob_interval,
  "\n",
  "Extreme-tail quantile/CDF round trip:",
  format(stable_cdf$p, digits = 16),
  "\n",
  sep = ""
)

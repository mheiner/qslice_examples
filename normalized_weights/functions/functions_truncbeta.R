## helper functions adapted from suggestions by ChatGPT Codex for numerical stability of decreasing weights example

.log_addexp <- function(lx, ly) {
  ## Calculates log(exp(lx) + exp(ly)).
  ## From R stats src/nmath/pgamma.c
  m <- pmax(lx, ly)
  out <- m + log1p(exp(-abs(lx - ly)))
  out[is.infinite(m) & m < 0.0] <- -Inf
  out
}

.log1mexp <- function(x) {
  ## Calculates log(1 - exp(x)) for x <= 0.
  stopifnot(all(x <= 0.0))

  out <- numeric(length(x))
  use_log1p <- x < -0.6931471805599453 # -0.6931472 = -log(2)

  out[use_log1p] <- log1p(-exp(x[use_log1p]))
  out[!use_log1p] <- log(-expm1(x[!use_log1p]))

  out
}

.log_subexp <- function(lx, ly) {
  ## Calculates log(exp(lx) - exp(ly)).
  ## From R stats src/nmath/pgamma.c
  if (any(ly > lx)) {
    stop("`ly` must not exceed `lx`.")
  }
  lx + .log1mexp(ly - lx)
}

.bounded_log_probs <- function(lb, ub, fn_cdf_free) {
  ## Log probabilities to the left of, within, and to the right of a
  ## truncation interval.
  stopifnot(
    length(lb) == 1L,
    length(ub) == 1L,
    is.finite(lb),
    is.finite(ub),
    lb < ub,
    is.function(fn_cdf_free)
  )

  log_half <- -0.6931471805599453

  lp_lb <- fn_cdf_free(lb, lower.tail = TRUE, log.p = TRUE)
  lsp_lb <- fn_cdf_free(lb, lower.tail = FALSE, log.p = TRUE)
  lp_ub <- fn_cdf_free(ub, lower.tail = TRUE, log.p = TRUE)
  lsp_ub <- fn_cdf_free(ub, lower.tail = FALSE, log.p = TRUE)

  ## Store log(F(lb)) and log(1 - F(ub)) using the smaller tail directly.
  lp_left <- if (lp_lb <= log_half) lp_lb else .log1mexp(lsp_lb)
  lp_right <- if (lsp_ub <= log_half) lsp_ub else .log1mexp(lp_ub)

  if (lp_ub <= log_half) {
    ## The interval is wholly in the lower-probability half.
    lp_interval <- .log_subexp(lp_ub, lp_left)
  } else if (lsp_lb <= log_half) {
    ## The interval is wholly in the upper-probability half.
    lp_interval <- .log_subexp(lsp_lb, lp_right)
  } else {
    ## The interval contains the median: subtract the two excluded tails.
    lp_interval <- .log1mexp(.log_addexp(lp_left, lp_right))
  }

  if (!is.finite(lp_interval)) {
    stop(
      "The truncation interval has zero probability at machine precision.",
      call. = FALSE
    )
  }

  list(
    lp_left = lp_left,
    lp_interval = lp_interval,
    lp_right = lp_right
  )
}

.beta_truncation_setup <- function(shape1, shape2, lower, upper) {
  stopifnot(
    length(shape1) == 1L,
    length(shape2) == 1L,
    length(lower) == 1L,
    length(upper) == 1L,
    is.finite(shape1),
    is.finite(shape2),
    shape1 > 0.0,
    shape2 > 0.0,
    !is.na(lower),
    !is.na(upper)
  )

  lower <- max(lower, 0.0)
  upper <- min(upper, 1.0)
  if (lower >= upper) {
    stop(
      "The truncation interval must have positive width within (0, 1).",
      call. = FALSE
    )
  }

  fn_cdf_free <- function(x, lower.tail = TRUE, log.p = FALSE) {
    pbeta(
      x,
      shape1 = shape1,
      shape2 = shape2,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  log_probs <- .bounded_log_probs(
    lb = lower,
    ub = upper,
    fn_cdf_free = fn_cdf_free
  )

  c(
    list(
      shape1 = shape1,
      shape2 = shape2,
      lower = lower,
      upper = upper
    ),
    log_probs
  )
}

qtruncbeta <- function(p, shape1, shape2, lower, upper) {
  setup <- .beta_truncation_setup(shape1, shape2, lower, upper)

  log_half <- -0.6931471805599453

  q_one <- function(p_now) {
    if (is.na(p_now) || p_now < 0.0 || p_now > 1.0) {
      stop("Probabilities in `p` must be in [0, 1].", call. = FALSE)
    }
    if (p_now == 0.0) {
      return(setup$lower)
    }
    if (p_now == 1.0) {
      return(setup$upper)
    }

    lp <- log(p_now)
    lsp <- log1p(-p_now)
    lp_free <- .log_addexp(
      setup$lp_left,
      lp + setup$lp_interval
    )

    if (lp_free <= log_half) {
      qbeta(
        lp_free,
        shape1 = setup$shape1,
        shape2 = setup$shape2,
        lower.tail = TRUE,
        log.p = TRUE
      )
    } else {
      lsp_free <- .log_addexp(
        setup$lp_right,
        lsp + setup$lp_interval
      )
      qbeta(
        lsp_free,
        shape1 = setup$shape1,
        shape2 = setup$shape2,
        lower.tail = FALSE,
        log.p = TRUE
      )
    }
  }

  q_out <- vapply(p, q_one, numeric(1))

  list(
    q = q_out,
    # Retained for compatibility/diagnostics only. This can legitimately
    # underflow to zero; never take log(prob_interval).
    prob_interval = exp(setup$lp_interval),
    # Use this value in every log-density calculation.
    log_prob_interval = setup$lp_interval
  )
}

ptruncbeta <- function(q, shape1, shape2, lower, upper) {
  setup <- .beta_truncation_setup(shape1, shape2, lower, upper)

  log_half <- -0.6931471805599453

  p_one <- function(q_now) {
    if (is.na(q_now)) {
      return(q_now)
    }
    if (q_now <= setup$lower) {
      return(0.0)
    }
    if (q_now >= setup$upper) {
      return(1.0)
    }

    lp_q <- pbeta(
      q_now,
      shape1 = setup$shape1,
      shape2 = setup$shape2,
      lower.tail = TRUE,
      log.p = TRUE
    )

    if (lp_q <= log_half) {
      lp_numerator <- .log_subexp(lp_q, setup$lp_left)
      out <- exp(lp_numerator - setup$lp_interval)
    } else {
      lsp_q <- pbeta(
        q_now,
        shape1 = setup$shape1,
        shape2 = setup$shape2,
        lower.tail = FALSE,
        log.p = TRUE
      )
      lp_upper <- .log_subexp(lsp_q, setup$lp_right) -
        setup$lp_interval
      out <- -expm1(lp_upper)
    }

    min(1.0, max(0.0, out))
  }

  list(
    p = vapply(q, p_one, numeric(1)),
    prob_interval = exp(setup$lp_interval),
    log_prob_interval = setup$lp_interval
  )
}

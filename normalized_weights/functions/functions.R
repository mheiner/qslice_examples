logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

log1m_exp <- function(lx) {
  if (lx > 0.0) {
    stop(paste0("log1m_exp input ", lx, "\nlx must be less than 0"))
  } else if (lx > -.Machine$double.eps) {
    out <- -Inf
    warning("log1m_exp input ", lx, "\nreturning -Inf")
  } else if (lx > -0.6931471805599453) {
    # .693147 ~= log(2.0)

    out <- log(-expm1(lx))
  } else {
    out <- log1p(-exp(lx))
  }

  out
}


w_to_v <- function(w, tol = 1e-10) {
  # weights to stick-breaking variables; returns same length as argument
  stopifnot(all(w >= 0.0), abs(sum(w) - 1.0) < tol)
  w / c(1.0, 1.0 - cumsum(w)[-length(w)])
}

v_to_w <- function(v) {
  exp(v_to_logw(v))
}

lw_to_v <- function(lw) {
  # log weights to stick-breaking variables; returns same length as argument
  K <- length(lw)

  if (K == 1L) {
    return(1.0)
  }

  # log_remaining[k] = log(sum(exp(lw[k:K])))
  log_remaining <- numeric(K)
  log_remaining[K] <- lw[K]

  for (k in seq.int(K - 1L, 1L)) {
    a <- lw[k]
    b <- log_remaining[k + 1L]

    # Stable scalar log(exp(a) + exp(b)).
    if (a >= b) {
      log_remaining[k] <- a + log1p(exp(b - a))
    } else {
      log_remaining[k] <- b + log1p(exp(a - b))
    }
  }

  lv <- lw - log_remaining

  # Mathematically lv <= 0. pmin protects against tiny positive
  # values caused by floating-point rounding.
  v <- exp(pmin(lv, 0.0))

  # Under the same-length convention, the terminal stick is exactly one.
  v[K] <- 1.0

  v
}

v_to_lw <- function(v) {
  # stick-breaking variables to weights; returns same length as argument
  stopifnot(
    is.numeric(v),
    length(v) >= 1L,
    all(is.finite(v)),
    all(v > 0.0),
    all(v <= 1.0),
    all(v[-length(v)] < 1.0)
  )

  m <- length(v)
  lv <- log(v)

  if (m == 1L) {
    return(lv)
  }

  # Only nonterminal complements enter the stick products.
  l1mv <- log1p(-v[-m])

  lv + c(0.0, cumsum(l1mv))
}

lw_to_logits <- function(lw) {
  # returns length one LESS than argument
  K <- length(lw)
  lw[1:(K - 1)] - lw[K]
}

w_to_logits <- function(w) {
  # returns length one LESS than argument
  stopifnot(all(w <= 1.0), all(w >= 0.0), sum(w) == 1.0)
  lw <- log(w)
  lw_to_logits(lw)
}

logits_to_lw <- function(logits) {
  # returns length one MORE than argument
  lnum <- c(logits, 0.0)
  ldenom <- logsumexp(lnum)
  lnum - ldenom
}

logits_to_w <- function(logits) {
  # returns length one MORE than argument
  exp(logits_to_lw(logits))
}

revcumsum <- function(x) rev(cumsum(rev(x[-1])))


## original target
log_targ_theta <- function(theta) {
  thetaC <- rep(theta, each = n) * CC
  denoms <- rowSums(thetaC)
  nums <- thetaC[cbind(1:n, x)]
  llik <- log(nums) - log(denoms)
  lpri <- (a0 - 1) * log(theta)

  sum(llik) + sum(lpri)
}

log_targ_theta_log <- function(ltheta) {
  ## also requires logC
  ## this is the original target density on theta computed using ltheta = log(theta)
  lthetalC <- rep(ltheta, each = n) + lCC
  ldenoms <- apply(lthetalC, 1L, logsumexp)
  lnums <- lthetalC[cbind(1:n, x)]
  llik <- lnums - ldenoms
  lpri <- (a0 - 1.0) * ltheta

  sum(llik) + sum(lpri)
}

## marginal target for logits = logit(theta)
log_targ_logits <- function(logits) {
  ltheta <- logits_to_lw(logits)
  out0 <- log_targ_theta_log(ltheta)
  ljacob <- sum(ltheta)
  out0 + ljacob
}

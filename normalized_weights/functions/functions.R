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
  } else if (lx > -0.693147) {
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
  # stick-breaking variables to weights; returns same length as argument
  stopifnot(all(v >= 0.0), all(v <= 1.0))

  ## below same as   # c(1.0, cumprod(1.0 - v[-length(v)])) * v
  lv <- log(v)
  l1mv <- sapply(lv, log1m_exp)
  lw <- c(0.0, cumsum(l1mv[-length(v)])) + lv
  exp(lw)
}

w_to_logits <- function(w) {
  # returns length one LESS than argument
  K <- length(w)
  log(w[1:(K - 1)]) - log(w[K])
}

logits_to_w <- function(logits) {
  # returns length one MORE than argument
  elogits <- exp(logits)
  sumelogits <- sum(elogits)
  sump1 <- sumelogits + 1.0
  c(elogits, 1.0) / sump1
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

## marginal target for logits = logit(theta)
log_targ_logits <- function(logits) {
  theta <- logits_to_w(logits)
  out0 <- log_targ_theta(theta)
  ljacob <- sum(log(theta))
  out0 + ljacob
}

slice_elliptical_mv_spherical <- function(z, sig, log_target) {
  nEvaluations <- 0
  k <- length(z)
  f <- function(z) {
    nEvaluations <<- nEvaluations + 1
    log_target(z)
  }
  fz <- f(z)$value
  stopifnot(fz > -Inf)
  y <- log(runif(1)) + fz
  nu <- rnorm(k, 0, sd = sig)
  twopi <- 2 * pi
  theta <- runif(1, 0, twopi)
  theta_min <- theta - twopi
  theta_max <- theta
  repeat {
    z1 <- z * cos(theta) + nu * sin(theta)
    fz1 <- f(z1)
    if (y < fz1$value) {
      return(list(z = z1, x = fz1$x, nEvaluations = nEvaluations))
    }
    if (theta < 0) {
      theta_min <- theta
    } else {
      theta_max <- theta
    }
    theta <- runif(1, theta_min, theta_max)
  }
}


slice_genelliptical_mv_logits <- function(
  x = NULL, # logits
  z = NULL, # transformed
  log_target_logits,
  mu,
  Sig,
  df,
  is_chol = FALSE,
  static_pseudo = FALSE, # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  type = "unordered"
) {
  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

  if (isTRUE(static_pseudo)) {
    # was given z
    k <- length(z)
  } else {
    # was given x
    k <- length(x)
    z <- drop(forwardsolve(SigL, (x - mu)))
  }

  stopifnot(length(mu) == k)
  stopifnot(dim(Sig) == c(k, k))

  a <- 0.5 * (df + k)
  b <- 0.5 * (df + sum(z^2))
  s <- 1.0 / rgamma(1, shape = a, rate = b)

  lff <- function(zz) {
    xx <- drop(SigL %*% zz + mu)
    ssz <- sum(zz^2)

    if (type == "unordered") {
      out <- log_target_logits(xx) + a * log1p(ssz / df)
    } else if (type == "decreasing") {
      if (all(diff(xx) < 0.0) & all(xx > 0.0)) {
        # xx refers to logits here
        out <- log_target_logits(xx) + a * log1p(ssz / df)
      } else {
        out <- -Inf
      }
    }

    list(value = out, z = zz, x = xx)
  }

  out_ess <- slice_elliptical_mv_spherical(
    z = z,
    sig = sqrt(s),
    log_target = lff
  )

  out_ess
}


mcmc_gess_logit <- function(
  n_iter,
  n_thin,
  state,
  log_target,
  mu,
  Sig,
  df,
  is_chol = FALSE,
  static_pseudo, # T/R: could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  type = "unordered",
  verbose = TRUE
) {
  Km1 <- length(state$logits)

  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

  sims_logits <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(sims_logits) <- paste0("logit", 1:Km1)
  simsZ <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(simsZ) <- paste0("z", 1:Km1)

  neval <- numeric(n_iter)

  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }

  if (isTRUE(static_pseudo)) {
    z_now <- forwardsolve(SigL, (state$logits - mu)) |> drop()
  } else {
    z_now <- NULL
  }

  for (i in 1:n_iter) {
    for (ii in 1:n_thin) {
      tmp <- slice_genelliptical_mv_logits(
        x = state$logits,
        z = z_now,
        log_target = log_target,
        mu = mu,
        Sig = SigL,
        df = df,
        is_chol = TRUE,
        static_pseudo = static_pseudo, # T/R: could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
        type = type
      )

      z_now <- tmp$z |> drop()
      state$logits <- tmp$x |> drop()
    }

    neval[i] <- tmp$nEvaluations
    sims_logits[i, ] <- state$logits
    simsZ[i, ] <- z_now

    if (isTRUE(verbose)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (isTRUE(verbose)) {
    close(pb)
  }

  list(state = state, sims_logits = sims_logits, simsZ = simsZ, neval = neval)
}

runif_equator <- function(theta) {
  ## draw from uniform on great subsphere on S_{d-2} with respect to x on S_{d-1} (Habeck, 2023)
  d <- length(theta)
  v1 <- rnorm(d)
  v2 <- v1 - (crossprod(theta, v1) %*% theta)
  v3 <- v2 / sqrt(sum(v2^2))
  drop(v3)
}

geodesic_shrinkage <- function(
  log_target,
  log_thresh,
  r_now,
  theta_now,
  SigL,
  mu
) {
  ## Implemented from Schar et al, 2023
  twopi <- 2.0 * pi
  y <- runif_equator(theta_now)
  omega_max <- runif(1, min = 0.0, max = twopi)
  omega_min <- omega_max - twopi

  accpt <- FALSE
  iter <- 1

  while (isFALSE(accpt) && iter < 5000) {
    omega_now <- runif(1, min = omega_min, max = omega_max)
    theta_out <- theta_now * cos(omega_now) + y * sin(omega_now)
    z_out <- drop(r_now * theta_out)
    x_out <- drop((SigL %*% z_out) + mu)
    if (log_target(x_out) > log_thresh) {
      accpt <- TRUE
    } else {
      # shrink
      if (omega_now < 0.0) {
        omega_min <- omega_now
      } else {
        omega_max <- omega_now
      }
      iter <- iter + 1
    }
  }
  list(x = theta_out, nEvaluations = iter)
}

stepsrhink_rz <- function(
  log_target,
  log_thresh,
  r_now,
  theta_now,
  d,
  SigL,
  mu,
  w_step
) {
  ## Based on Schar et al 2023
  ## log_target here will refer only to g(x(z)) while z(r, theta) = r*theta

  n_eval_now <- 0
  ltarg <- function(x) {
    n_eval_now <<- n_eval_now + 1 # so we can count while stepping out
    log_target(x)
  }

  U <- runif(1, min = 0.0, max = 1.0)
  r_min <- max(r_now - U * w_step, 0.0)
  r_max <- r_now + (1.0 - U) * w_step

  z_rmin <- drop(r_min * theta_now)
  x_rmin <- (SigL %*% z_rmin) + mu

  ## step
  while ((r_min > 0) && (ltarg(x_rmin) + (d - 1) * log(r_min) > log_thresh)) {
    r_min <- max(r_min - w_step, 0.0)
    z_rmin <- drop(r_min * theta_now)
    x_rmin <- (SigL %*% z_rmin) + mu
  }

  z_rmax <- drop(r_max * theta_now)
  x_rmax <- (SigL %*% z_rmax) + mu

  while ((ltarg(x_rmax) + (d - 1) * log(r_max)) > log_thresh) {
    r_max <- r_max + w_step
    z_rmax <- drop(r_max * theta_now)
    x_rmax <- (SigL %*% z_rmax) + mu
  }

  r_cand <- runif(1, min = r_min, max = r_max)
  z_cand <- drop(r_cand * theta_now)
  x_cand <- drop((SigL %*% z_cand) + mu)

  ## shrink
  while ((ltarg(x_cand) + (d - 1) * log(r_cand)) <= log_thresh) {
    if (r_cand < r_now) {
      r_min <- r_cand
    } else {
      r_max <- r_cand
    }
    r_cand <- runif(1, min = r_min, max = r_max)
    z_cand <- drop(r_cand * theta_now)
    x_cand <- drop((SigL %*% z_cand) + mu)
  }

  list(x = x_cand, z = z_cand, r = r_cand, nEvaluations = n_eval_now)
}

slice_tgpolar_logits <- function(
  log_target, # evaluates on original x (logits)
  logits = NULL,
  thetaz = NULL,
  rz = NULL,
  static_pseudo = FALSE, # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  mu,
  SigL,
  w_step_rz,
  type = "unordered"
) {
  ## Gibbsian Polar Slice Sampling (Schar, Habeck, Rudolf, 2023) with location-scale transformation
  ## Adapted from Algorithm 1 on Schar, Habeck, Rudolph (2023)
  ## log_target evaluated on original x (logits)

  if (isTRUE(static_pseudo)) {
    d <- length(thetaz)
    z_old <- drop(thetaz * rz)
    x_old <- (SigL %*% z_old) + mu
  } else {
    d <- length(logits)
    x_old <- logits
    z_old <- forwardsolve(SigL, (x_old - mu)) |> drop()
    rz <- sqrt(sum(z_old^2))
    thetaz <- z_old / rz
  }
  dlrz_old <- (d - 1) * log(rz)

  lthresh <- log(runif(1, min = 0.0, max = 1.0)) + log_target(x_old) + dlrz_old

  if (type == "decreasing") {
    log_target_use <- function(logits) {
      if (all(diff(logits) < 0.0) & all(logits > 0.0)) {
        out <- log_target(logits)
      } else {
        out <- -Inf
      }
    }
  } else if (type == "unordered") {
    log_target_use <- log_target # rename function
  }

  tmp <- geodesic_shrinkage(
    log_target = log_target_use,
    log_thresh = lthresh - dlrz_old, # offset this one because the function only compares target to threshold
    r_now = rz,
    theta_now = thetaz,
    SigL = SigL,
    mu = mu
  )

  neval_theta <- tmp$nEvaluations
  thetaz_new <- tmp$x

  tmp <- stepsrhink_rz(
    log_target = log_target_use,
    log_thresh = lthresh, # don't offset this one because this one compares the full f1 to threshold
    r_now = rz,
    theta_now = thetaz_new,
    d = d,
    SigL = SigL,
    mu = mu,
    w_step = w_step_rz
  )

  neval_r <- tmp$nEvaluations
  rz_new <- tmp$r
  z_new <- tmp$z
  x_new <- tmp$x # stepshrink_rz built to be second update to save on one calculation of z and x

  list(
    thetaz = thetaz_new,
    rz = rz_new,
    z = z_new,
    x = x_new,
    nEvaluations = 1 + neval_theta + neval_r,
    nEvaluations_theta = neval_theta,
    nEvaluations_r = neval_r
  )
}


mcmc_gpss_logit <- function(
  n_iter,
  n_thin,
  state,
  log_target, # evaluated as x (logit)
  mu,
  Sig,
  is_chol = FALSE,
  static_pseudo = FALSE, # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  w_step_rz = 1.0,
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

  neval <- numeric(n_iter)

  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }

  if (isTRUE(static_pseudo)) {
    z_now <- forwardsolve(SigL, (state$logits - mu)) |> drop() # can just run chain on Z without having to recompute Z
    rz_now <- sqrt(sum(z_now^2))
    thetaz_now <- z_now / rz_now
  } else {
    thetaz_now <- NULL # will be calculated at each iter with new pseudo
    rz_now <- NULL # will be calculated at each iter with new pseudo
  }

  for (i in 1:n_iter) {
    for (ii in 1:n_thin) {
      tmp <- slice_tgpolar_logits(
        log_target = log_target,
        logits = state$logits,
        thetaz = thetaz_now,
        rz = rz_now,
        static_pseudo = static_pseudo,
        mu = mu,
        SigL = SigL,
        w_step_rz = w_step_rz,
        type = type
      )

      state$logits <- drop(tmp$x)
      thetaz_now <- drop(tmp$thetaz)
      rz_now <- tmp$rz
    }

    neval[i] <- tmp$nEvaluations
    sims_logits[i, ] <- state$logits
    simsZ[i, ] <- tmp$z

    if (isTRUE(verbose)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (isTRUE(verbose)) {
    close(pb)
  }

  list(state = state, sims_logits = sims_logits, simsZ = simsZ, neval = neval)
}

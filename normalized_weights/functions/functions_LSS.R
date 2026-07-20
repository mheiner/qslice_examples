slice_latent_logit <- function(
  logits,
  z, # z is forwardsolve(SigL, logits - mu)
  log_target, # evaluates logit
  mu,
  SigL,
  static_pseudo = static_pseudo, # T/F: could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  latent_s, # latent parameters, same length as z
  latent_scale, # parameter for auxiliary distribution of latent_s
  type = "unordered"
) {
  if (type == "decreasing") {
    log_target_use <- function(logits) {
      if (all(diff(logits) < 0.0) & all(logits > 0.0)) {
        out <- log_target(logits)
      } else {
        out <- -Inf
      }
    }
  } else if (type == "unordered") {
    log_target_use <- log_target # copy function
  }

  if (isTRUE(static_pseudo)) {
    k <- length(z)
    logits <- drop(SigL %*% z + mu)
  } else {
    k <- length(logits)
    z <- forwardsolve(SigL, logits - mu) |> drop()
  }

  lfx <- log_target_use(logits)
  nEvaluations <- 1

  if (is.nan(lfx) | is.infinite(lfx) | is.na(lfx)) {
    stop(paste0(
      "input z = ",
      paste(z, collapse = " "),
      "\n",
      "log_target(x(z)) = ",
      lfx,
      "\n"
    ))
  }

  y <- log(runif(1, min = 0.0, max = 1.0)) + lfx

  half_s <- latent_s / 2.0
  ll <- runif(k, min = z - half_s, max = z + half_s)
  s_now <- -latent_scale * log(runif(k, min = 0.0, max = 1.0)) + 2 * abs(ll - z)
  half_s_now <- s_now / 2.0

  L <- ll - half_s_now
  R <- ll + half_s_now

  repeat {
    z1 <- L + runif(k, min = 0.0, max = 1.0) * (R - L)
    logits1 <- drop(SigL %*% z1 + mu)

    lfx1 <- log_target_use(logits1)
    nEvaluations <- nEvaluations + 1

    if (is.nan(lfx1) | is.na(lfx1)) {
      warning(paste0(
        "input z = ",
        paste(z, collapse = " "),
        "\n",
        "proposed z1 = ",
        paste(z1, collapse = " "),
        "\n",
        "log_target(x(z1)) = ",
        lfx1,
        "\n"
      ))

      lfx1
    }

    if (y < lfx1) {
      return(list(
        z = z1,
        logits = logits1,
        latent_s = s_now,
        nEvaluations = nEvaluations
      ))
    }

    shrinkL <- (z1 < z)

    if (any(shrinkL)) {
      indx <- which(shrinkL)
      L[indx] <- z1[indx]
    }

    if (any(!shrinkL)) {
      indx <- which(!shrinkL)
      R[indx] <- z1[indx]
    }
  }
}


mcmc_lss_logit <- function(
  n_iter,
  n_thin,
  state, # state includes latent_s vector
  log_target, # input log-target function evaluated on logit = SigL %*% z + mu scale
  mu,
  Sig,
  is_chol = FALSE,
  static_pseudo = FALSE, # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  latent_scale = 1.0, # z standardizes the variance of each dimension, so we go with a single tuning parameter
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
  simsS <- matrix(NA, ncol = Km1, nrow = n_iter)

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
      tmp <- slice_latent_logit(
        logits = state$logits,
        z = z_now, # z is forwardsolve(SigL, logits - mu)
        log_target = log_target, # evaluates logit
        mu = mu,
        SigL = SigL,
        static_pseudo = static_pseudo, # T/F: could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
        latent_s = state$latent_s, # latent parameters, same length as z
        latent_scale = latent_scale, # parameter for auxiliary distribution of latent_s
        type = type
      )

      z_now <- tmp$z
      state$logits <- tmp$logits
      state$latent_s <- tmp$latent_s
    }

    neval[i] <- tmp$nEvaluations
    sims_logits[i, ] <- state$logits
    simsZ[i, ] <- z_now
    simsS[i, ] <- state$latent_s

    if (isTRUE(verbose)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (isTRUE(verbose)) {
    close(pb)
  }

  list(
    state = state,
    sims_logits = sims_logits,
    simsZ = simsZ,
    simsS = simsS,
    neval = neval
  )
}

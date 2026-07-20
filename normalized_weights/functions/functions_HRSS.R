slice_hyperrect_logit <- function(
  logits = NULL,
  z = NULL, # forwardsolve(SigL, logits - mu)
  log_target, # evaluates logit
  mu,
  SigL,
  static_pseudo = FALSE, # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  w = NULL,
  L = NULL,
  R = NULL,
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

  if (!all(is.null(w))) {
    L <- z - runif(k) * w
    R <- L + w
  } else if (all(is.null(c(L, R)))) {
    L <- rep(0.0, k)
    R <- rep(1.0, k)
  }

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
      return(list(z = z1, logits = logits1, nEvaluations = nEvaluations))
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


mcmc_hrss_logit <- function(
  n_iter,
  n_thin,
  state,
  log_target, # input log-target function evaluated on logit = SigL %*% z + mu scale
  mu,
  Sig,
  is_chol = FALSE,
  static_pseudo = FALSE, # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  w = 2.0, # z standardizes the variance of each dimension, so we go with a single tuning parameter
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
    z_now <- forwardsolve(SigL, (state$logits - mu)) |> drop()
  } else {
    z_now <- NULL
  }

  for (i in 1:n_iter) {
    for (ii in 1:n_thin) {
      tmp <- slice_hyperrect_logit(
        logits = state$logits,
        z = z_now,
        log_target = log_target,
        mu = mu,
        SigL = SigL,
        static_pseudo = static_pseudo, # T/F could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
        w = w,
        L = NULL,
        R = NULL,
        type = type
      )

      z_now <- tmp$z
      state$logits <- tmp$logits
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

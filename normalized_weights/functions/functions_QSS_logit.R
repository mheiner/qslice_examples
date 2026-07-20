slice_hyperrect_indep_t <- function(
  u, # takes u
  log_target_u # evaluates on u
) {
  k <- length(u)

  lhu <- log_target_u(u)$value
  nEvaluations <- 1

  if (is.infinite(lhu) | is.nan(lhu) | is.na(lhu)) {
    stop(paste0(
      "input u = ",
      paste(x, collapse = " "),
      "\n",
      "log_target_u(u) = ",
      lhu,
      "\n"
    ))
  }

  y <- log(runif(1, min = 0.0, max = 1.0)) + lhu

  L <- rep(0.0, k)
  R <- rep(1.0, k)

  repeat {
    u1 <- L + runif(k, min = 0.0, max = 1.0) * (R - L)

    tmp1 <- log_target_u(u1)
    lhu1 <- tmp1$value
    nEvaluations <- nEvaluations + 1

    if (is.nan(lhu1) | is.na(lhu1)) {
      stop(paste0(
        "input u = ",
        paste(u, collapse = " "),
        "\n",
        "proposed u1 = ",
        paste(u1, collapse = " "),
        "\n",
        "proposed z = ",
        paste(tmp1$z, collapse = " "),
        "\n",
        "proposed logits = ",
        paste(tmp1$x, collapse = " "),
        "\n",
        "log_target_u(u1) = ",
        lhu1,
        "\n"
      ))
    }

    if (nEvaluations > 3001) {
      stop(paste0(
        "Exceeded maximum shrinkage steps for U\n",
        "current L: ",
        L,
        "\n",
        "current R: ",
        R,
        "\n",
        "current proposed u1: ",
        u1,
        "\n"
      ))
    }

    if (y < lhu1) {
      return(list(u = u1, x = tmp1$x, z = tmp1$z, nEvaluations = nEvaluations))
    }

    shrinkL <- (u1 < u)

    if (any(shrinkL)) {
      indx <- which(shrinkL)
      L[indx] <- u1[indx]
    }

    if (any(!shrinkL)) {
      indx <- which(!shrinkL)
      R[indx] <- u1[indx]
    }
  }
}


slice_hyperrect_logits <- function(
  x = NULL, # logits
  z = NULL, # transformed
  u = NULL, # CDF-transformed
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
    # work with incoming z and u
    k <- length(z)
  } else {
    # work with incoming x (logits), need to recompute z and u
    k <- length(x)
    z <- drop(forwardsolve(SigL, (x - mu)))
    u <- pt(z, df = df)
  }

  stopifnot(length(mu) == k)
  stopifnot(dim(Sig) == c(k, k))

  a <- 0.5 * (df + k)
  b <- 0.5 * (df + sum(z^2))

  lff <- function(uu) {
    zz <- qt(uu, df = df)
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

    list(value = out, x = xx, z = zz, u = uu)
  }

  out_qss <- slice_hyperrect_indep_t(u = u, log_target_u = lff)

  out_qss
}


mcmc_qss_logit <- function(
  n_iter,
  n_thin,
  state,
  log_target, # log_target_logits
  mu,
  Sig,
  df,
  is_chol = FALSE,
  static_pseudo, # T/F: could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  type = "unordered",
  verbose = TRUE
) {
  ## Pseudo-target is assumed independent student ts AFTER rotational transformation.

  ## With QSS, transformation to U is always necessary. However, if pseudo-target is
  ## not static, U must be recomputed from the current z and pseudo.

  Km1 <- length(state$logits)

  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

  simsU <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(simsU) <- paste0("U", 1:Km1)

  simsZ <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(simsZ) <- paste0("z", 1:Km1)

  sims_logits <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(sims_logits) <- paste0("logit", 1:Km1)

  neval <- numeric(n_iter)

  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }

  z_now <- drop(forwardsolve(SigL, (state$logits - mu)))
  u_now <- pt(z_now, df = df)

  for (i in 1:n_iter) {
    for (ii in 1:n_thin) {
      tmp <- slice_hyperrect_logits(
        x = state$logits,
        z = z_now,
        u = u_now,
        log_target = log_target,
        mu = mu,
        Sig = SigL,
        df = df,
        is_chol = TRUE,
        static_pseudo = static_pseudo, # T/F: could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
        type = type
      )

      u_now <- tmp$u |> drop()
      z_now <- tmp$z |> drop()
      state$logits <- tmp$x |> drop()
    }

    neval[i] <- tmp$nEvaluations
    sims_logits[i, ] <- state$logits
    simsZ[i, ] <- z_now
    simsU[i, ] <- u_now

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
    simsU = simsU,
    neval = neval
  )
}

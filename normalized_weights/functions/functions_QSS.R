.decseq_upper <- function(v_previous) {
  # min(v / (1 - v), 1), avoiding overflow and unnecessary division.
  if (v_previous >= 0.5) 1.0 else v_previous / (1.0 - v_previous)
}

V_decseq_from_U <- function(u_vec, shape1_vec, shape2_vec) {
  Km1 <- length(u_vec)
  K <- Km1 + 1L

  stopifnot(
    Km1 >= 1L,
    length(shape1_vec) == Km1,
    length(shape2_vec) == Km1,
    all(is.finite(shape1_vec)),
    all(is.finite(shape2_vec)),
    all(shape1_vec > 0.0),
    all(shape2_vec > 0.0),
    all(!is.na(u_vec)),
    all(u_vec >= 0.0),
    all(u_vec <= 1.0)
  )

  lower <- 1.0 / (K - seq_len(Km1) + 1.0)
  upper <- numeric(Km1)
  vv <- numeric(Km1)
  log_prob_interval <- numeric(Km1)

  upper[1] <- 1.0
  for (k in seq_len(Km1)) {
    if (k > 1L) {
      upper[k] <- .decseq_upper(vv[k - 1L])
    }
    if (lower[k] >= upper[k]) {
      # warning(
      #   "Empty decreasing-stick interval at component ",
      #   k,
      #   ": lower = ",
      #   lower[k],
      #   ", upper = ",
      #   upper[k],
      #   "."
      # )
      return(list(
        valid = FALSE,
        failure_component = k,
        lower = lower,
        upper = upper
      ))
    }

    tmp <- qtruncbeta(
      p = u_vec[k],
      shape1 = shape1_vec[k],
      shape2 = shape2_vec[k],
      lower = lower[k],
      upper = upper[k]
    )
    vv[k] <- tmp$q
    log_prob_interval[k] <- tmp$log_prob_interval
  }

  list(
    v = vv,
    lv = log(vv),
    lower = lower,
    upper = upper,
    log_prob_intervals = log_prob_interval,
    valid = TRUE
  )
}

U_from_Vdecseq <- function(v_vec, shape1_vec, shape2_vec) {
  Km1 <- length(v_vec)
  K <- Km1 + 1L

  stopifnot(
    Km1 >= 1L,
    length(shape1_vec) == Km1,
    length(shape2_vec) == Km1,
    all(is.finite(shape1_vec)),
    all(is.finite(shape2_vec)),
    all(shape1_vec > 0.0),
    all(shape2_vec > 0.0),
    all(is.finite(v_vec))
  )

  lower <- 1.0 / (K - seq_len(Km1) + 1.0)
  upper <- c(
    1.0,
    vapply(v_vec[-Km1], .decseq_upper, numeric(1))
  )
  u <- numeric(Km1)
  log_prob_interval <- numeric(Km1)

  for (k in seq_len(Km1)) {
    if (v_vec[k] <= lower[k] || v_vec[k] >= upper[k]) {
      stop(
        "Stick component ",
        k,
        " is outside its open conditional support: v = ",
        v_vec[k],
        ", interval = (",
        lower[k],
        ", ",
        upper[k],
        ").",
        call. = FALSE
      )
    }

    tmp <- ptruncbeta(
      q = v_vec[k],
      shape1 = shape1_vec[k],
      shape2 = shape2_vec[k],
      lower = lower[k],
      upper = upper[k]
    )
    u[k] <- tmp$p
    log_prob_interval[k] <- tmp$log_prob_interval
  }

  list(
    u = u,
    lower = lower,
    upper = upper,
    log_prob_intervals = log_prob_interval
  )
}


slice_hyperrect_stick <- function(
  u, # takes u
  log_target_u # evaluates on u and returns eval and v
) {
  k <- length(u)

  lhu <- log_target_u(u)$lhu
  nEvaluations <- 1L

  if (is.infinite(lhu) | is.nan(lhu) | is.na(lhu)) {
    stop(paste0(
      "input u = ",
      paste(u, collapse = " "),
      "\n",
      "log_target_u(u) = ",
      lhu,
      "\n"
    ))
  }

  y <- log(runif(1, min = 0.0, max = 1.0)) + lhu

  L <- rep(0.0, k)
  R <- rep(1.0, k)

  nAttempts <- 0L

  repeat {
    nAttempts <- nAttempts + 1L
    u1 <- L + runif(k, min = 0.0, max = 1.0) * (R - L)
    tmp1 <- log_target_u(u1)
    lhu1 <- tmp1$lhu
    if (isTRUE(tmp1$eval)) {
      nEvaluations <- nEvaluations + 1L
    }

    if (is.nan(lhu1) | is.na(lhu1)) {
      stop(paste0(
        "input u = ",
        paste(u, collapse = " "),
        "\n",
        "proposed u1 = ",
        paste(u1, collapse = " "),
        "\n",
        "log_target_u(u1) = ",
        lhu1,
        "\n"
      ))
    }

    if (nAttempts > 3001) {
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
      return(list(u = u1, v = tmp1$v, nEvaluations = nEvaluations))
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


slice_quantile_stick <- function(
  v,
  u = NULL,
  log_target_theta_log,
  pseudo,
  static_pseudo = FALSE,
  type = "unordered"
) {
  ## Pseudo-target is assumed beta. pseudo is a list with the shape params

  if (isFALSE(static_pseudo)) {
    # we need to recompute U every iter
    if (type == "unordered") {
      u <- pbeta(v, shape1 = pseudo$shape1, shape2 = pseudo$shape2)
    } else if (type == "decreasing") {
      UfromV <- U_from_Vdecseq(
        v,
        shape1_vec = pseudo$shape1,
        shape2_vec = pseudo$shape2
      )
      u <- UfromV$u
    }
  }

  if (type == "unordered") {
    lhu <- function(u) {
      vv <- qbeta(u, shape1 = pseudo$shape1, shape2 = pseudo$shape2)

      Km1 <- length(u)

      lvv <- log(vv)
      l1mv <- log1p(-vv)
      cs_l1mv <- cumsum(l1mv)

      if (any(is.infinite(c(lvv, l1mv)))) {
        # out of support: return -Inf

        warning(paste0(
          "QSS lhu input u = ",
          paste(u, collapse = " "),
          "\n",
          "log(v) = ",
          paste(lvv, collapse = " "),
          "\n",
          "log(1-v) = ",
          paste(l1mv, collapse = " "),
          "\n",
          "Out of support. Proceeding by returning -Inf"
        ))
        out <- -Inf
        eval <- FALSE
      } else {
        ltheta <- c(0.0, cs_l1mv[-Km1]) + lvv # length K - 1
        ltheta_K <- cs_l1mv[Km1]

        ltheta <- c(ltheta, ltheta_K) # becomes length K

        lt0 <- log_target_theta_log(ltheta)
        eval <- TRUE

        ## straightforward incorporation of Jacobian
        # ljacob <- sum(ltheta) - sum(lvv)
        # lt <- lt0 + ljacob
        #
        # ldens_pseudo <- sum((pseudo$shape1 - 1.0)*lvv) + sum((pseudo$shape2 - 1.0)*l1mv)
        # out <- lt - ldens_pseudo

        ## tricky incorporation of Jacobian
        ldens_pseudo <- sum((pseudo$shape1 - 1.0) * lvv) +
          sum((pseudo$shape2 - Km1:1) * l1mv) # includes Jacobian from numerator
        out <- lt0 - ldens_pseudo # dividing by this pseudo puts Jacobian to numerator
      }

      list(lhu = out, v = vv, eval = eval)
    }
  } else if (type == "decreasing") {
    lhu <- function(u) {
      Km1 <- length(u)
      VfromU <- V_decseq_from_U(
        u_vec = u,
        shape1_vec = pseudo$shape1,
        shape2_vec = pseudo$shape2
      )

      if (isTRUE(VfromU$valid)) {
        vv <- VfromU$v
        lvv <- VfromU$lv
        l1mv <- log1p(-vv)
        cs_l1mv <- cumsum(l1mv)

        if (any(is.infinite(c(lvv, l1mv)))) {
          # out of support: return -Inf

          warning(paste0(
            "QSS lhu input u = ",
            paste(u, collapse = " "),
            "\n",
            "log(v) = ",
            paste(lvv, collapse = " "),
            "\n",
            "log(1-v) = ",
            paste(l1mv, collapse = " "),
            "\n",
            "Out of support. Proceeding by returning -Inf"
          ))

          out <- -Inf
          eval <- FALSE
        } else {
          ltheta <- c(0.0, cs_l1mv[-Km1]) + lvv # length K - 1
          ltheta_K <- cs_l1mv[Km1]

          ltheta <- c(ltheta, ltheta_K) # becomes length K

          lt0 <- log_target_theta_log(ltheta)
          eval <- TRUE

          ## straightforward incorporation of Jacobian theta to V
          # ljacob <- sum(ltheta) - sum(lvv)
          # lt <- lt0 + ljacob
          #
          # ldens_pseudo <- sum((pseudo$shape1 - 1.0)*lvv) + sum((pseudo$shape2 - 1.0)*l1mv) - sum(VfromU$log_prob_intervals[-1])
          # out <- lt - ldens_pseudo

          ## tricky incorporation of Jacobian theta to V
          ldens_pseudo <- sum((pseudo$shape1 - 1.0) * lvv) +
            sum((pseudo$shape2 - Km1:1) * l1mv) -
            sum(VfromU$log_prob_intervals[-1]) # includes Jacobian
          out <- lt0 - ldens_pseudo # dividing by this pseudo puts Jacobian to numerator
        }
      } else {
        out <- -Inf
        vv <- NA
        eval <- FALSE
      }
      list(lhu = out, v = vv, eval = eval)
    }
  }

  shyp <- slice_hyperrect_stick(u, log_target_u = lhu)

  list(u = shyp$u, v = shyp$v, nEvaluations = shyp$nEvaluations)
}


mcmc_qss_stick <- function(
  n_iter,
  n_thin,
  state,
  log_target_theta_log,
  pseudo,
  static_pseudo = TRUE,
  type = "unordered",
  verbose = TRUE
) {
  ## Pseudo-target is assumed beta. pseudo is a list with the shape params.

  ## With QSS, transformation to U is always necessary. However, if pseudo-target is
  ## not static, U must be recomputed from the current V and pseudo.

  Km1 <- length(state$v)

  simsU <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(simsU) <- paste0("U", 1:Km1)

  simsV <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(simsV) <- paste0("v", 1:Km1)

  neval <- numeric(n_iter)

  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }

  if (type == "unordered") {
    u_now <- pbeta(state$v, shape1 = pseudo$shape1, shape2 = pseudo$shape2)
  } else if (type == "decreasing") {
    UfromV <- U_from_Vdecseq(
      state$v,
      shape1_vec = pseudo$shape1,
      shape2_vec = pseudo$shape2
    )
    u_now <- UfromV$u
  }

  for (i in 1:n_iter) {
    for (ii in 1:n_thin) {
      tmp <- slice_quantile_stick(
        v = state$v,
        u = u_now,
        log_target_theta_log = log_target_theta_log,
        pseudo = pseudo,
        static_pseudo = static_pseudo,
        type = type
      )
      u_now <- tmp$u
      state$v <- tmp$v
    }

    neval[i] <- tmp$nEvaluations
    simsU[i, ] <- u_now
    simsV[i, ] <- state$v

    if (isTRUE(verbose)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (isTRUE(verbose)) {
    close(pb)
  }

  list(state = state, simsV = simsV, simsU = simsU, neval = neval)
}

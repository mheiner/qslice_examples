qtruncbeta <- function(p, shape1, shape2, lower, upper) {
  p_lower <- pbeta(lower, shape1 = shape1, shape2 = shape2)
  p_upper <- pbeta(upper, shape1 = shape1, shape2 = shape2)
  prob_interval <- p_upper - p_lower
  p_out <- p * prob_interval + p_lower
  q_out <- qbeta(p_out, shape1 = shape1, shape2 = shape2)
  list(q = q_out, prob_interval = prob_interval)
}

V_decseq_from_U <- function(u_vec, shape1_vec, shape2_vec) {
  Km1 <- length(u_vec)
  K <- Km1 + 1
  lower <- numeric(Km1)
  upper <- numeric(Km1)
  vv <- numeric(Km1)
  lvv <- numeric(Km1)
  prob_interval <- numeric(Km1)

  lower[1] <- 1.0 / K
  upper[1] <- 1.0

  p_lower1 <- pbeta(lower[1], shape1 = shape1_vec[1], shape2 = shape2_vec[1])
  prob_interval[1] <- 1.0 - p_lower1

  vv[1] <- qbeta(
    u_vec[1] * prob_interval[1] + p_lower1,
    shape1 = shape1_vec[1],
    shape2 = shape2_vec[1]
  )
  lvv[1] <- log(vv[1])

  for (k in 2:Km1) {
    lower[k] <- 1.0 / (K - k + 1)
    upper[k] <- exp(lvv[k - 1] - log1m_exp(lvv[k - 1]))
    tmp <- qtruncbeta(
      p = u_vec[k],
      shape1 = shape1_vec[k],
      shape2 = shape2_vec[k],
      lower = lower[k],
      upper = upper[k]
    )
    vv[k] <- tmp$q
    lvv[k] <- log(tmp$q)
    prob_interval[k] <- tmp$prob_interval
  }

  if (any(is.nan(vv)) | any(is.na(vv))) {
    warning(paste0(
      "Decreasing V sequence returning Na or NaN\n",
      "V = ",
      paste(vv, collapse = " "),
      "\n"
    ))
  }

  if (any(is.nan(log(prob_interval))) | any(is.na(prob_interval))) {
    warning(paste0(
      "Na or NaN in probability intervals for truncated betas.\n",
      "V = ",
      paste(vv, collapse = " "),
      "\n"
    ))
  }

  list(
    v = vv,
    lv = lvv,
    lower = lower,
    upper = upper,
    prob_intervals = prob_interval
  )
}

U_from_Vdecseq <- function(v_vec, shape1_vec, shape2_vec) {
  Km1 <- length(v_vec)
  K <- Km1 + 1
  lower <- 1.0 / (K - 1:Km1 + 1)
  upper <- c(1.0, v_vec[-Km1] / (1.0 - v_vec[-Km1]))
  prob_interval <- pbeta(upper, shape1 = shape1_vec, shape2 = shape2_vec) -
    pbeta(lower, shape1 = shape1_vec, shape2 = shape2_vec)

  p0 <- pbeta(v_vec, shape1 = shape1_vec, shape2 = shape2_vec)
  p_lower <- pbeta(lower, shape1 = shape1_vec, shape2 = shape2_vec)

  u <- (p0 - p_lower) / prob_interval

  list(u = u, lower = lower, upper = upper, prob_intervals = prob_interval)
}


slice_hyperrect_stick <- function(
  u, # takes u
  log_target_u # evaluates on u and returns eval and v
) {
  k <- length(u)

  lhu <- log_target_u(u)$lhu
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
    lhu1 <- tmp1$lhu
    nEvaluations <- nEvaluations + 1

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
  log_target_theta,
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
      l1mv <- sapply(lvv, log1m_exp)

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
      } else {
        ltheta <- c(0.0, cumsum(l1mv[-Km1])) + lvv # will be K - 1
        ltheta_K <- sum(l1mv)

        theta <- exp(c(ltheta, ltheta_K))

        lt0 <- log_target_theta(theta)

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

      list(lhu = out, v = vv)
    }
  } else if (type == "decreasing") {
    lhu <- function(u) {
      Km1 <- length(u)
      VfromU <- V_decseq_from_U(
        u_vec = u,
        shape1_vec = pseudo$shape1,
        shape2 = pseudo$shape2
      )

      vv <- VfromU$v
      lvv <- VfromU$lv
      l1mv <- sapply(lvv, log1m_exp)

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
      } else {
        ltheta <- c(0.0, cumsum(l1mv[-Km1])) + lvv # will be K - 1
        ltheta_K <- sum(l1mv)

        theta <- exp(c(ltheta, ltheta_K))

        lt0 <- log_target_theta(theta)

        ## straightforward incorporation of Jacobian theta to V
        # ljacob <- sum(ltheta) - sum(lvv)
        # lt <- lt0 + ljacob
        #
        # ldens_pseudo <- sum((pseudo$shape1 - 1.0)*lvv) + sum((pseudo$shape2 - 1.0)*l1mv) - sum(log(VfromU$prob_intervals[-1]))
        # out <- lt - ldens_pseudo

        ## tricky incorporation of Jacobian theta to V
        ldens_pseudo <- sum((pseudo$shape1 - 1.0) * lvv) +
          sum((pseudo$shape2 - Km1:1) * l1mv) -
          sum(log(VfromU$prob_intervals[-1])) # includes Jacobian
        out <- lt0 - ldens_pseudo # dividing by this pseudo puts Jacobian to numerator
      }

      list(lhu = out, v = vv)
    }
  }

  shyp <- slice_hyperrect_stick(u, log_target_u = lhu)

  list(u = shyp$u, v = shyp$v, nEvaluations = shyp$nEvaluations)
}


mcmc_qss_stick <- function(
  n_iter,
  n_thin,
  state,
  log_target_theta,
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
        log_target_theta = log_target_theta,
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

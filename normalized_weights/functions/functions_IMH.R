imh_beta <- function(vv, log_target_theta, pseudo) {
  Km1 <- length(vv)

  lvv <- log(vv)
  l1mv <- sapply(lvv, log1m_exp)

  ltheta <- c(0.0, cumsum(l1mv[-Km1])) + lvv # will be K - 1
  ltheta_K <- sum(l1mv)

  theta <- exp(c(ltheta, ltheta_K))

  lt0 <- log_target_theta(theta)

  ## straightforward incorporation of Jacobian to theta (remains commented out)
  # ljacob <- sum(ltheta) - sum(lvv)
  # lt <- lt0 + ljacob
  # ldens_pseudo <- sum((pseudo$shape1 - 1.0)*lvv) + sum((pseudo$shape2 - 1.0)*l1mv)
  # out <- lt - ldens_pseudo

  ## tricky incorporation of Jacobian to theta
  ldens_pseudo <- sum((pseudo$shape1 - 1.0) * lvv) +
    sum((pseudo$shape2 - Km1:1) * l1mv) # includes Jacobian
  lfx_old <- lt0 - ldens_pseudo # dividing by this pseudo puts Jacobian to numerator

  vv_cand <- rbeta(Km1, shape1 = pseudo$shape1, shape2 = pseudo$shape2)

  lvv_cand <- log(vv_cand)
  l1mv_cand <- sapply(lvv_cand, log1m_exp)

  if (any(is.infinite(c(lvv_cand, l1mv_cand)))) {
    # out of support: return -Inf

    warning(paste0(
      "IMH proposal\n",
      "log(v) = ",
      paste(lvv_cand, collapse = " "),
      "\n",
      "log(1-v) = ",
      paste(l1mv_cand, collapse = " "),
      "\n",
      "Out of support. Proceeding by returning -Inf"
    ))

    lfx_cand <- -Inf
  } else {
    ltheta_cand <- c(0.0, cumsum(l1mv_cand[-Km1])) + lvv_cand # will be K - 1
    ltheta_K_cand <- sum(l1mv_cand)

    theta_cand <- exp(c(ltheta_cand, ltheta_K_cand))

    lt0_cand <- log_target_theta(theta_cand)
    ldens_pseudo_cand <- sum((pseudo$shape1 - 1.0) * lvv_cand) +
      sum((pseudo$shape2 - Km1:1) * l1mv_cand) # includes Jacobian
    lfx_cand <- lt0_cand - ldens_pseudo_cand
  }

  lprob_accpt <- lfx_cand - lfx_old
  lu_accpt <- log(runif(1, min = 0.0, max = 1.0))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- vv_cand
  } else {
    out <- vv
  }
  out
}


imh_decbetaU <- function(
  v = NULL,
  u = NULL,
  log_target_theta,
  pseudo,
  static_pseudo = FALSE # could pseudo change from iteration to iteration? If not, we can build the chain on U and save computation
) {
  if (isTRUE(static_pseudo)) {
    Km1 <- length(u)
    tnx <- V_decseq_from_U(
      u_vec = u,
      shape1_vec = pseudo$shape1,
      shape2 = pseudo$shape2
    )
    v <- tnx$v
    lv <- tnx$lv
  } else {
    Km1 <- length(v)
    lv <- log(v)
    tnx <- U_from_Vdecseq(
      v_vec = v,
      shape1_vec = pseudo$shape1,
      shape2_vec = pseudo$shape2
    )
    u <- tnx$u
  }

  l1mv <- sapply(lv, log1m_exp)

  ltheta <- c(0.0, cumsum(l1mv[-Km1])) + lv # will be K - 1
  ltheta_K <- sum(l1mv)

  theta <- exp(c(ltheta, ltheta_K))

  lt0 <- log_target_theta(theta)

  ## tricky incorporation of Jacobian to theta
  ldens_pseudo <- sum((pseudo$shape1 - 1.0) * lv) +
    sum((pseudo$shape2 - Km1:1) * l1mv) -
    sum(log(tnx$prob_intervals[-1])) # includes Jacobian v to theta and prob_intervals of proposal
  lfx_old <- lt0 - ldens_pseudo # dividing by this pseudo puts Jacobian to numerator

  u_cand <- runif(Km1, min = 0.0, max = 1.0)
  VfromU_cand <- V_decseq_from_U(
    u_vec = u_cand,
    shape1_vec = pseudo$shape1,
    shape2 = pseudo$shape2
  )

  v_cand <- VfromU_cand$v
  lv_cand <- VfromU_cand$lv
  l1mv_cand <- sapply(lv_cand, log1m_exp)

  if (any(is.infinite(c(lv_cand, l1mv_cand)))) {
    # out of support: return -Inf

    warning(paste0(
      "IMH proposal\n",
      "log(v) = ",
      paste(lv_cand, collapse = " "),
      "\n",
      "log(1-v) = ",
      paste(l1mv_cand, collapse = " "),
      "\n",
      "Out of support. Proceeding by returning -Inf"
    ))

    lfx_cand <- -Inf
  } else {
    ltheta_cand <- c(0.0, cumsum(l1mv_cand[-Km1])) + lv_cand # will be K - 1
    ltheta_K_cand <- sum(l1mv_cand)

    theta_cand <- exp(c(ltheta_cand, ltheta_K_cand))

    lt0_cand <- log_target_theta(theta_cand)
    ldens_pseudo_cand <- sum((pseudo$shape1 - 1.0) * lv_cand) +
      sum((pseudo$shape2 - Km1:1) * l1mv_cand) -
      sum(log(VfromU_cand$prob_intervals[-1])) # includes Jacobian
    lfx_cand <- lt0_cand - ldens_pseudo_cand
  }

  lprob_accpt <- lfx_cand - lfx_old
  lu_accpt <- log(runif(1, min = 0.0, max = 1.0))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- list(u = u_cand, v = v_cand, theta = theta_cand)
  } else {
    out <- list(u = u, v = v, theta = theta)
  }
  out
}


mcmc_imh_stick <- function(
  n_iter,
  n_thin,
  state,
  log_target_theta,
  pseudo,
  static_pseudo = FALSE,
  type = "unordered",
  verbose = TRUE
) {
  Km1 <- length(state$v)

  simsV <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(simsV) <- paste0("v", 1:Km1)

  neval <- rep(2, n_iter)

  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }

  if (type == "decreasing") {
    # working in U is mandatory
    UfromV <- U_from_Vdecseq(
      v_vec = state$v,
      shape1_vec = pseudo$shape1,
      shape2_vec = pseudo$shape2
    )
    u_now <- UfromV$u

    simsU <- matrix(NA, ncol = Km1, nrow = n_iter)
    colnames(simsU) <- paste0("u", 1:Km1)
  } else {
    simsU <- NULL
  }

  for (i in 1:n_iter) {
    for (ii in 1:n_thin) {
      if (type == "unordered") {
        state$v <- imh_beta(
          state$v,
          log_target_theta = log_target_theta,
          pseudo = pseudo
        ) # no need for transformation
      } else if (type == "decreasing") {
        # working in U is mandatory
        tmp <- imh_decbetaU(
          v = state$v,
          u = u_now,
          log_target_theta = log_target_theta,
          pseudo = pseudo,
          static_pseudo = static_pseudo # T/F: could pseudo change from iteration to iteration? If not, we can build the chain on U and save computation
        )
        u_now <- tmp$u
        state$v <- tmp$v
      }
    }

    simsV[i, ] <- state$v
    if (type == "decreasing") {
      simsU[i, ] <- u_now
    }

    if (isTRUE(verbose)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (isTRUE(verbose)) {
    close(pb)
  }

  list(state = state, simsV = simsV, simsU = simsU, neval = neval)
}

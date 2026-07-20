imh_logits <- function(
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
  ## Since this evaluates on Z, we cannot enforce "decreasing" constraint in the proposal.
  ## If we did enforce it, we would need to sample logits (actually, not even that would work),
  ## which would require recomputing the proposal density (and inverting sigL)
  ## at each iteration.

  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

  if (isTRUE(static_pseudo)) {
    # work with incoming z
    k <- length(z)
  } else {
    # work with incoming x (logits), need to recompute z
    k <- length(x)
    z <- drop(forwardsolve(SigL, (x - mu)))
  }

  stopifnot(length(mu) == k)
  stopifnot(dim(Sig) == c(k, k))

  a <- 0.5 * (df + k)
  b <- 0.5 * (df + sum(z^2))

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

    list(value = out, x = xx, z = zz)
  }

  z_cand <- rt(k, df = df)

  lrat_old <- lff(z)
  lrat_cand <- lff(z_cand)

  lprob_accpt <- lrat_cand$value - lrat_old$value
  lu_accpt <- log(runif(1, min = 0.0, max = 1.0))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- lrat_cand
  } else {
    out <- lrat_old
  }

  out
}

mcmc_imh_logit <- function(
  n_iter,
  n_thin,
  state,
  log_target_logits,
  mu,
  Sig,
  df,
  is_chol = FALSE,
  static_pseudo = FALSE,
  type = "unordered",
  verbose = TRUE
) {
  Km1 <- length(state$logits)

  if (isTRUE(is_chol)) {
    SigL <- Sig
  } else {
    SigL <- t(chol(Sig))
  }

  simsZ <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(simsZ) <- paste0("z", 1:Km1)

  sims_logits <- matrix(NA, ncol = Km1, nrow = n_iter)
  colnames(sims_logits) <- paste0("logit", 1:Km1)

  neval <- rep(2, n_iter)

  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }

  z_now <- drop(forwardsolve(SigL, (state$logits - mu)))

  for (i in 1:n_iter) {
    for (ii in 1:n_thin) {
      tmp <- imh_logits(
        x = state$logits,
        z = z_now,
        log_target_logits = log_target_logits,
        mu = mu,
        Sig = SigL,
        df = df,
        is_chol = TRUE,
        static_pseudo = static_pseudo, # T/F? could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
        type = type
      )

      z_now <- tmp$z |> drop()
      state$logits <- tmp$x |> drop()
    }

    simsZ[i, ] <- z_now
    sims_logits[i, ] <- state$logits

    if (isTRUE(verbose)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (isTRUE(verbose)) {
    close(pb)
  }

  list(state = state, simsZ = simsZ, sims_logits = sims_logits, neval = neval)
}

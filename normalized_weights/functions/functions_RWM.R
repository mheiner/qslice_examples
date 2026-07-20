rw_logit <- function(
  logits,
  z, # forwardsolve(SigL, logits - mu)
  log_target, # evaluates logit
  mu,
  SigL,
  static_pseudo = FALSE, # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  Cscale, # tuning param
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

  z1 <- z + Cscale * rnorm(k)
  logits1 <- drop(SigL %*% z1 + mu)
  lfx1 <- log_target_use(logits1)

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

  lprob_accpt <- lfx1 - lfx
  lu_accpt <- log(runif(1, min = 0.0, max = 1.0))

  if (isTRUE(lu_accpt < lprob_accpt)) {
    out <- list(z = z1, logits = logits1, nEvaluations = 2, accpt = TRUE)
  } else {
    out <- list(z = z, logits = logits, nEvaluations = 2, accpt = FALSE)
  }

  out
}


mcmc_rw_logit <- function(
  n_iter,
  n_thin,
  state,
  log_target, # input log-target function evaluated on logit = SigL %*% z + mu scale
  mu,
  Sig,
  is_chol = FALSE,
  static_pseudo = FALSE, # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
  Cscale = 1.0, # z standardizes the variance of each dimension, so we go with a single tuning parameter
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
  naccpt <- 0

  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }

  if (isTRUE(static_pseudo)) {
    z_now <- forwardsolve(SigL, (state$logits - mu)) |> drop()
    logits_now <- state$logits
  } else {
    z_now <- NULL
    logits_now <- state$logits
  }

  for (i in 1:n_iter) {
    for (ii in 1:n_thin) {
      tmp <- rw_logit(
        logits = logits_now,
        z = z_now, # z is forwardsolve(SigL, logits - mu)
        log_target = log_target, # evaluates logit
        mu = mu,
        SigL = SigL,
        static_pseudo = static_pseudo, # T/F: could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
        Cscale = Cscale, # tuning param
        type = type
      )

      z_now <- tmp$z
      logits_now <- tmp$logits
      state$logits <- tmp$logits
    }

    sims_logits[i, ] <- logits_now
    simsZ[i, ] <- z_now
    naccpt <- naccpt + tmp$accpt
    neval[i] <- tmp$nEvaluations

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
    neval = neval,
    naccpt = naccpt
  )
}

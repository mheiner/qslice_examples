mcmc_sample <- function(
  n_iter,
  n_thin,
  state,
  sampler,
  log_target_theta,
  log_target_logits,
  pseudo,
  target_type,
  tuning_params,
  verbose = FALSE
) {
  if (sampler == "QSS_stick") {
    tt <- system.time(
      mc <- mcmc_qss_stick(
        n_iter = n_iter,
        n_thin = n_thin,
        state = state,
        log_target_theta = log_targ_theta,
        pseudo = list(shape1 = pseudo$shape1, shape2 = pseudo$shape2),
        static_pseudo = pseudo$static,
        type = target_type,
        verbose = verbose
      )
    )
  } else if (sampler == "QSS_logit") {
    tt <- system.time(
      mc <- mcmc_qss_logit(
        n_iter = n_iter,
        n_thin = n_thin,
        state = state,
        log_target = log_targ_logits, # log_target_logits
        mu = pseudo$mu,
        Sig = pseudo$Sig,
        df = pseudo$df,
        is_chol = FALSE,
        static_pseudo = pseudo$static,
        type = target_type,
        verbose = verbose
      )
    )
  } else if (sampler == "IMH_stick") {
    tt <- system.time(
      mc <- mcmc_imh_stick(
        n_iter = n_iter,
        n_thin = n_thin,
        state = state,
        log_target_theta = log_targ_theta,
        pseudo = list(shape1 = pseudo$shape1, shape2 = pseudo$shape2),
        static_pseudo = pseudo$static,
        type = target_type,
        verbose = verbose
      )
    )
  } else if (sampler == "IMH_logit") {
    tt <- system.time(
      mc <- mcmc_imh_logit(
        n_iter = n_iter,
        n_thin = n_thin,
        state = state,
        log_target_logits = log_targ_logits,
        mu = pseudo$mu,
        Sig = pseudo$Sig,
        df = pseudo$df,
        is_chol = FALSE,
        static_pseudo = pseudo$static,
        type = target_type,
        verbose = verbose
      )
    )
  } else if (sampler == "GESS") {
    tt <- system.time(
      mc <- mcmc_gess_logit(
        n_iter = n_iter,
        n_thin = n_thin,
        state = state,
        log_target = log_targ_logits,
        mu = pseudo$mu,
        Sig = pseudo$Sig,
        df = pseudo$df,
        is_chol = FALSE,
        static_pseudo = pseudo$static,
        type = target_type,
        verbose = verbose
      )
    )
  } else if (sampler == "GPSS") {
    tt <- system.time(
      mc <- mcmc_gpss_logit(
        n_iter = n_iter,
        n_thin = n_thin,
        state = state,
        log_target = log_targ_logits,
        mu = pseudo$mu,
        Sig = pseudo$Sig,
        is_chol = FALSE,
        static_pseudo = pseudo$static,
        w_step_rz = tuning_params[["w_step_rz"]], # z standardizes the variance of each dimension, so we go with a single tuning parameter
        type = target_type,
        verbose = verbose
      )
    )
  } else if (sampler == "HRSS") {
    tt <- system.time(
      mc <- mcmc_hrss_logit(
        n_iter = n_iter,
        n_thin = n_thin,
        state = state,
        log_target = log_targ_logits,
        mu = pseudo$mu,
        Sig = pseudo$Sig,
        is_chol = FALSE,
        static_pseudo = pseudo$static,
        w = tuning_params[["w"]], # z standardizes the variance of each dimension, so we go with a single tuning parameter
        type = target_type,
        verbose = verbose
      )
    )
  } else if (sampler == "LSS") {
    tt <- system.time(
      mc <- mcmc_lss_logit(
        n_iter = n_iter,
        n_thin = n_thin,
        state = state,
        log_target = log_targ_logits,
        mu = pseudo$mu,
        Sig = pseudo$Sig,
        is_chol = FALSE,
        static_pseudo = pseudo$static,
        latent_scale = tuning_params[["latent_scale"]], # z standardizes the variance of each dimension, so we go with a single tuning parameter
        type = target_type,
        verbose = verbose
      )
    )
  } else if (sampler == "RWM") {
    tt <- system.time(
      mc <- mcmc_rw_logit(
        n_iter = n_iter,
        n_thin = n_thin,
        state = state,
        log_target = log_targ_logits,
        mu = pseudo$mu,
        Sig = pseudo$Sig,
        is_chol = FALSE,
        static_pseudo = pseudo$static,
        Cscale = tuning_params[["Cscale"]], # z standardizes the variance of each dimension, so we go with a single tuning parameter
        type = target_type,
        verbose = verbose
      )
    )
  } else {
    stop("No sampler selected.")
  }

  if (sampler %in% c("QSS_stick", "IMH_stick")) {
    # sampled stick-breaking v (sims outputs v, simsU included)

    SigL <- t(chol(pseudo$Sig))
    mc$sims_logits <- matrix(NA, nrow = n_iter, ncol = length(mc$state$logits))
    mc$simsZ <- matrix(NA, nrow = n_iter, ncol = length(mc$state$logits))

    for (ii in 1:n_iter) {
      theta_now <- v_to_w(c(mc$simsV[ii, ], 0.9999999999))
      mc$sims_logits[ii, ] <- w_to_logits(theta_now)
      mc$simsZ[ii, ] <- forwardsolve(
        SigL,
        (mc$sims_logits[ii, ] - pseudo$mu)
      ) |>
        drop()
    }

    mc$state$theta <- v_to_w(c(mc$state$v, 0.9999999999))
    mc$state$logits <- w_to_logits(mc$state$theta)
  } else {
    # sampled logits (sims outputs logits, simsZ included)

    mc$state$theta <- logits_to_w(mc$state$logits)
    mc$state$v <- w_to_v(mc$state$theta)[1:(length(mc$state$theta) - 1)]

    mc$simsV <- matrix(NA, nrow = n_iter, ncol = length(mc$state$v))

    for (ii in 1:n_iter) {
      theta_now <- logits_to_w(mc$sims_logits[ii, ])
      mc$simsV[ii, ] <- w_to_v(theta_now)[1:(length(theta_now) - 1)]
    }
  }

  list(mc = mc, time = tt)
}

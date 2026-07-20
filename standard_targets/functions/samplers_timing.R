############## Random Walk ################

random_walk_sampler <- function(state, n_iter, lf, support, c) {
  draws <- numeric(n_iter)
  n.accept <- 0

  for (i in 1:n_iter) {
    x_cand <- rnorm(1, mean = state$x, sd = c)

    if (x_cand >= support[1] && x_cand <= support[2]) {
      logr <- lf(x_cand) - lf(state$x)

      u <- runif(1, min = 0.0, max = 1.0)

      if (log(u) < logr) {
        state$x <- x_cand
        n.accept <- n.accept + 1
      }
    }

    draws[i] <- state$x
  }

  list(
    draws = draws,
    counter = 2 * n_iter,
    n.accept = n.accept,
    n_iter = n_iter,
    state = state
  )
}


############ Stepping Out Eval #############

## stepping out and shrinkage procedure of Neal (2003)

stepping_out_sampler <- function(state, n_iter, lf, w, max) {
  counter <- 0
  draws <- numeric(n_iter)

  for (i in 1:n_iter) {
    out <- qslice::slice_stepping_out(
      x = state$x,
      log_target = lf,
      w = w,
      max = max
    )
    state$x <- out$x
    draws[i] <- out$x
    counter <- counter + out$nEvaluations
  }

  list(draws = draws, counter = counter, n_iter = n_iter, state = state)
}


################ GESS ##############

## generalized elliptical slice sampler (Nishihara, 2014)

gess_sampler <- function(state, n_iter, lf, mu, sigma, degf) {
  counter <- 0
  draws <- numeric(n_iter)

  for (i in 1:n_iter) {
    out <- qslice::slice_genelliptical(
      x = state$x,
      log_target = lf,
      mu = mu,
      sigma = sigma,
      df = degf
    )

    state$x <- out$x
    draws[i] <- out$x
    counter <- counter + out$nEvaluations
  }

  list(draws = draws, counter = counter, n_iter = n_iter, state = state)
}


############ Latent Eval ############

# latent slice sampler (Li and Walker, 2023)

latent_sampler <- function(state, n_iter, lf, rate) {
  counter <- 0
  draws <- latent_s <- numeric(n_iter)

  for (i in 1:n_iter) {
    out <- qslice::slice_latent(
      x = state$x,
      s = state$s,
      log_target = lf,
      rate = rate
    )
    state$x <- out$x
    state$s <- out$s
    draws[i] <- out$x
    latent_s[i] <- out$s
    counter <- counter + out$nEvaluations
  }

  list(
    draws = draws,
    latent_s = latent_s,
    counter = counter,
    n_iter = n_iter,
    state = state
  )
}


############## Quantile Slice Eval ################

quantile_sampler <- function(state, n_iter, lf, pseudo) {
  counter <- 0
  draws <- numeric(n_iter)
  Udraws <- numeric(n_iter)

  for (i in 1:n_iter) {
    out <- qslice::slice_quantile(x = state$x, log_target = lf, pseudo = pseudo)
    state$x <- out$x
    draws[i] <- out$x
    Udraws[i] <- out$u
    counter <- counter + out$nEvaluations
  }

  list(
    draws = draws,
    Udraws = Udraws,
    counter = counter,
    n_iter = n_iter,
    state = state
  )
}


############## Independence Metropolis Hastings ################

IMH_sampler <- function(state, n_iter, lf, pseudo) {
  draws <- numeric(n_iter)
  n.accept <- 0

  for (i in 1:n_iter) {
    tmp <- imh_pseudo(x = state$x, log_target = lf, pseudo = pseudo)
    state$x <- tmp$x
    draws[i] <- tmp$x
    n.accept <- tmp$accpt
  }

  list(
    draws = draws,
    counter = 2 * n_iter,
    n.accept = n.accept,
    n_iter = n_iter,
    state = state
  )
}


### universal timer

# function to evaluate Random Walk
sampler_time_eval <- function(
  type,
  state,
  n_iter,
  lf_func,
  support,
  settings,
  ess_log
) {
  if (type == "rw") {
    time <- system.time({
      mcmc_out <- random_walk_sampler(
        state = state,
        n_iter = n_iter,
        lf = lf_func,
        support = support,
        c = settings$tune_param
      )
    })
  } else if (type == "stepping") {
    time <- system.time({
      mcmc_out <- stepping_out_sampler(
        state = state,
        n_iter = n_iter,
        lf = lf_func,
        w = settings$tune_param,
        max = Inf
      )
    })
  } else if (type == "gess") {
    time <- system.time({
      mcmc_out <- gess_sampler(
        state = state,
        n_iter = n_iter,
        lf = lf_func,
        mu = settings$loc,
        sigma = settings$sc,
        degf = settings$degf
      )
    })
  } else if (type == "latent") {
    time <- system.time({
      mcmc_out <- latent_sampler(
        state = state,
        n_iter = n_iter,
        lf = lf_func,
        rate = settings$tune_param
      )
    })
  } else if (type == "Qslice") {
    time <- system.time({
      mcmc_out <- quantile_sampler(
        state = state,
        n_iter = n_iter,
        lf = lf_func,
        pseudo = settings
      )
    })
  } else if (type == "imh") {
    time <- system.time({
      mcmc_out <- IMH_sampler(
        state = state,
        n_iter = n_iter,
        lf = lf_func,
        pseudo = settings
      )
    })
  }

  if (isTRUE(ess_log)) {
    draws_ess <- log(mcmc_out$draws)
  } else {
    draws_ess <- mcmc_out$draws
  }

  ESS <- min(n_iter, coda::effectiveSize(coda::as.mcmc(draws_ess))) # not necessary

  out <- list()
  out$tbl <- data.frame(
    nEval = mcmc_out$counter,
    EffSamp = ESS,
    userTime = time['user.self'],
    sysTime = time['sys.self'],
    elapsedTime = time['elapsed']
  )

  out$draws <- mcmc_out$draws
  out$state <- mcmc_out$state

  out
}

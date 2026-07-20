args <- commandArgs(trailingOnly = TRUE)
ii <- as.numeric(args[1]) # job id
dte <- as.numeric(args[2])

##### for testing
# ii <- 1650
# dte <- 260704
#####

library("coda")

source("functions/functions.R")
source("functions/functions_QSS.R")
source("functions/functions_QSS_logit.R")
source("functions/functions_IMH.R")
source("functions/functions_IMH_logit.R")
source("functions/functions_GESS.R")
source("functions/functions_GPSS.R")
source("functions/functions_HRSS.R")
source("functions/functions_LSS.R")
source("functions/functions_RWM.R")
source("functions/sampler_super.R")
source("functions/tune.R")

sessionInfo()

load(paste0("input/schedule_all_", dte, ".rda"))
# head(sched, n = 15)

run_id <- job_order[ii]
(run_info <- sched[run_id, ])
(type <- as.character(run_info$target))
(sampler <- as.character(run_info$sampler))
(pseudo_type <- as.character(run_info$pseudo))
(K <- run_info$K)
(rep_id <- run_info$rep)

a0 <- 0.5 # theta ~ Dirichlet(a0)
n <- 500 # sample size
C_sig <- 2.0 # C lognormal parameter
C_mu <- 0.0 # C lognormal parameter

cat("whether to tune?\n")
(to_tune <- sampler %in% c("GPSS", "HRSS", "LSS", "RWM"))

set.seed(rep_id) # rep_id identifies each data set
source("0_sim_data.R")
cat("data:\n")
theta
cbind(xtab, theta)

n_burn <- 10e3
n_iter <- 50e3
n_chains <- 3

set.seed(dte + ii)


## set up data-based pseudo-target
cat("Setting up data-based pseudo-target\n")
ensure_ergod <- TRUE

(ps_naive <- a0 + xtab)
(theta_hat_naive <- ps_naive / sum(ps_naive))
(expos <- colSums(CC / drop(CC %*% theta_hat_naive)))
(theta_hat_expos0 <- xtab / expos)
(theta_hat_expos <- theta_hat_expos0 / sum(theta_hat_expos0))
(ps_expos <- a0 + theta_hat_expos * n)

pseudo_scale <- 0.9

if (isTRUE(ensure_ergod)) {
  ps_expos_rat <- min(ps_naive / ps_expos)
  pseudo_ref <- ps_expos * ps_expos_rat * pseudo_scale

  cbind(ps_expos, ps_naive, pseudo_ref, diff = ps_naive - pseudo_ref)
  stopifnot(all((ps_naive - pseudo_ref) >= 0.0))
} else {
  pseudo_ref <- ps_expos * pseudo_scale
}

pseudo_ref_rcs <- revcumsum(pseudo_ref)

pseudo_ref_use <- pseudo_ref[1:(K - 1)]
pseudo_ref_rcs_use <- pseudo_ref_rcs[1:(K - 1)]

cat("Pseudo beta params:\n")
cbind(pseudo_ref_use, pseudo_ref_rcs_use)

## Aitchison 1986, p 127 K-L approx for pseudo on logits
cat("Pseudo logit:\n")
(cand_mu <- digamma(ps_expos[-K]) - digamma(ps_expos[K]))
(cand_Sig <- matrix(trigamma(ps_expos[K]), nrow = K - 1, ncol = K - 1) +
  diag(trigamma(ps_expos[-K])))

pseudo <- list(
  shape1 = pseudo_ref_use,
  shape2 = pseudo_ref_rcs_use,
  mu = cand_mu,
  Sig = cand_Sig,
  df = ifelse(sampler %in% c("QSS_logit", "IMH_logit"), 10, 5),
  static = TRUE # will change according to settings
)

tuning_params <- list(w_step_rz = 1.0, w = 2, latent_scale = 1.0, Cscale = 0.6)

## initialize one chain
init_state <- function() {
  state <- list()
  (state$v <- rbeta(K - 1, shape1 = pseudo$shape1, shape2 = pseudo$shape2))
  (state$theta <- v_to_w(c(state$v, 0.9999999999)))
  if (type == "decreasing") {
    state$theta <- sort(state$theta, decreasing = TRUE)
    state$v <- w_to_v(state$theta)[-K]
  }
  (state$logits <- w_to_logits(state$theta)) # must be all decreasing and greater than 0 if type == "decreasing"
  (state$latent_s <- rgamma(
    K - 1,
    shape = 2.0,
    scale = tuning_params$latent_scale
  )) # for Latent Slice
  state
}

state <- init_state()

if (isTRUE(to_tune) || pseudo_type == "samples") {
  ## get initial samples (initial samples MUST be of reasonable quality, so init samples should use GESS or QSS above)
  cat("Initial sampling\n")
  mc <- mcmc_sample(
    n_iter = ifelse(type == "decreasing", 10e3, 2e3),
    n_thin = 1,
    state = state,
    sampler = "GESS",
    log_target_theta = log_mtarg_theta,
    log_target_logits = log_mtarg_logits,
    pseudo = pseudo,
    target_type = type,
    tuning_params = tuning_params,
    verbose = FALSE
  )

  mc <- mcmc_sample(
    n_iter = 2000,
    n_thin = ifelse(type == "decreasing", 20, 4),
    state = mc$mc$state,
    sampler = "GESS",
    log_target_theta = log_mtarg_theta,
    log_target_logits = log_mtarg_logits,
    pseudo = pseudo,
    target_type = type,
    tuning_params = tuning_params,
    verbose = FALSE
  )
}

## update the pseudo-target
if (pseudo_type == "samples") {
  cat("Sample-based pseudo:\n")

  pseudo_scale_tuned <- 0.9

  if (sampler %in% c("QSS_stick", "IMH_stick")) {
    # tuning on sticks

    qvm <- colMeans(mc$mc$simsV)
    qvv <- apply(mc$mc$simsV, 2, var)
    qvM <- qvm * (1.0 - qvm) / qvv - 1.0

    (pseudo_ref_use <- qvm * qvM * pseudo_scale_tuned)
    (pseudo_ref_rcs_use <- ((1.0 - qvm) * qvM * pseudo_scale_tuned))

    pseudo$shape1 <- pseudo_ref_use
    pseudo$shape2 <- pseudo_ref_rcs_use

    ensure_ergod <- FALSE

    ## will also need mu and Sig on logits to compute Z for ESS

    (pseudo$mu <- colMeans(mc$mc$sims_logits))
    (pseudo$Sig <- cov(mc$mc$sims_logits) / pseudo_scale_tuned)
  } else {
    # tuning on logits

    (pseudo$mu <- colMeans(mc$mc$sims_logits))
    (pseudo$Sig <- cov(mc$mc$sims_logits) / pseudo_scale_tuned)
  }

  pseudo$static <- TRUE
} else if (pseudo_type == "data") {
  pseudo$static <- FALSE
}

## tune extra params in certain samplers
if (isTRUE(to_tune)) {
  cat("Tuning:\n")

  if (sampler == "RWM") {
    tune_param_name <- "Cscale"
  }
  if (sampler == "HRSS") {
    tune_param_name <- "w"
  }
  if (sampler == "GPSS") {
    tune_param_name <- "w_step_rz"
  }
  if (sampler == "LSS") {
    tune_param_name <- "latent_scale"
  }

  tune_bnds_init <- c(
    0.5 * tuning_params[[tune_param_name]],
    8.0 * tuning_params[[tune_param_name]]
  )

  mc_tune <- tune(
    state = mc$mc$state,
    log_target = log_mtarg_logits, # each of the tuned algorithms work on logit space
    target_type = type,
    sampler = sampler,
    settings = tuning_params, # named list with tuning parameter values
    bnds_init = tune_bnds_init,
    pseudo = pseudo,
    param_name = tune_param_name,
    n_iter = 1000,
    n_grid = 3,
    n_rep = 2,
    n_rounds = 5,
    range_frac = 0.8,
    verbose = TRUE
  )

  # plot(mc_tune$lvals, mc_tune$lesps)
  # abline(v = log(mc_tune$val_opt), lty = 2)

  tuning_params[[tune_param_name]] <- mc_tune$val_opt
} else {
  tune_param_name <- "other"
  tuning_params[[tune_param_name]] <- NA
}


## initialize chains
states <- lapply(1:n_chains, function(j) init_state())


## burn in all chains
mc <- lapply(1:n_chains, function(j) {
  mcmc_sample(
    n_iter = n_burn,
    n_thin = 1,
    state = states[[j]],
    sampler = sampler,
    log_target_theta = log_mtarg_theta,
    log_target_logits = log_mtarg_logits,
    pseudo = pseudo,
    target_type = type,
    tuning_params = tuning_params,
    verbose = FALSE
  )
})

# plot(as.mcmc(mc[[1]]$mc$simsZ[,c(1,5,9)]))

for (j in 1:n_chains) {
  states[[j]] <- mc[[j]]$mc$state
}

## timing runs
mc <- lapply(1:n_chains, function(j) {
  mcmc_sample(
    n_iter = n_iter,
    n_thin = 1,
    state = states[[j]],
    sampler = sampler,
    log_target_theta = log_mtarg_theta,
    log_target_logits = log_mtarg_logits,
    pseudo = pseudo,
    target_type = type,
    tuning_params = tuning_params,
    verbose = FALSE
  )
})

## calculate stats
cat("Calculating stats:\n")
(tt <- sapply(mc, function(x) x$time["user.self"]))
(ESSZ <- sapply(mc, function(x) coda::effectiveSize(as.mcmc(x$mc$simsZ))))
(ESSZ_means <- pmax(colMeans(ESSZ), 1.0))
sum(ESSZ_means) / sum(tt)
mean(ESSZ_means / tt)
(ESSZ_mins <- apply(ESSZ, 2, min))
(ESpSZ_min_avg <- mean(ESSZ_mins / tt))
(Rhat_u95_avg <- tryCatch(
  {
    coda::gelman.diag(
      lapply(mc, function(x) as.mcmc(x$mc$simsZ)),
      confidence = 0.95,
      transform = FALSE,
      autoburnin = FALSE
    )$psrf[, 2] |>
      mean()
  },
  error = function(e) {
    cat("Warning: cannot compute R-hat\n")
    NA
  },
  finally = function() {
    NA
  }
))
(mean_neval <- sapply(mc, function(x) mean(x$mc$neval)))
(espitZ_avg <- mean(ESSZ / n_iter))
(IATZ_avg <- mean(n_iter / ESSZ_means))

logits_true <- w_to_logits(theta)
RMSE_logits_avg <- sapply(mc, function(x) {
  sqrt(mean((x$mc$sims_logits - rep(logits_true, each = n_iter))^2))
}) |>
  mean()


## report
tempDf <- data.frame(
  target = type,
  K = K,
  sampler = sampler,
  pseudo = pseudo_type,
  rep = run_info$rep,
  run_id = run_id,
  ii = ii,
  unif_ergod = ensure_ergod,
  tune_param_name = tune_param_name,
  tune_param_val = tuning_params[[tune_param_name]],
  n_chains = n_chains,
  n_iter = n_iter,
  nEval_mean = mean(mean_neval),
  ESSZ_mean_total = sum(ESSZ_means),
  userTime_total = sum(tt),
  ESpSZ_min_avg = ESpSZ_min_avg,
  Rhat_u95_avg = Rhat_u95_avg,
  espitZ_avg = espitZ_avg,
  IATZ_avg = IATZ_avg,
  RMSE_logits_avg = RMSE_logits_avg
)

## save summary
write.table(
  tempDf,
  file = paste0(
    "output/target",
    type,
    "_K",
    K,
    "_sampler",
    sampler,
    "_pseudo",
    pseudo_type,
    "_rep",
    run_info$rep,
    "_dte",
    dte,
    ".csv"
  ),
  append = FALSE,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)

print('finished')

quit(save = "no")

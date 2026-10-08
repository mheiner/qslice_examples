args <- commandArgs(trailingOnly = TRUE)
ii <- as.numeric(args[1]) # job id
dte <- as.numeric(args[2])

##### for testing
# ii <- 1964
# dte <- 240520
#####

library("qslice")
library("coda")
source("0_data.R")
source("0_prior.R")
source("functions/MH_samplers.R")
source("functions/full_conditionals.R")
source("functions/mcmc_gprior.R")
source("functions/tune.R")
source("functions/Laplace_g.R")
source("functions/invgamma.R")

load("schedule_all.rda")

sessionInfo()

n_iter <- 30e3

run_id <- job_order[ii]
(run_info <- sched[run_id, ])
target <- as.character(run_info$target)
logG <- grepl("log", target)
type <- as.character(run_info$type)
subtype <- as.character(run_info$subtype)

ess_log <- TRUE # measure effective sample size on log(g) regardless of chosen transformation

state0 <- list(iter = 0)
state0$psi <- prior$a_0 / prior$b_0
state0$g <- prior$g_max / 3
state0$beta <- prior$beta_0
state0

sampler0 <- list()

sampler0$g <- list(
  type = "stepping",
  subtype = NA,
  logG = logG,
  w = ifelse(logG, 5.0, 50.0),
  doextra = FALSE
)
sampler_tuned <- sampler0 # make a copy

### initial burn-in
mc_out <- mcmc_gprior(
  state = state0,
  prior = prior,
  data = dat,
  sampler = sampler0,
  n_iter = 2e3,
  save = FALSE,
  prog = 0
)

(state1 <- mc_out$state)

mc_out <- mcmc_gprior(
  state = state1,
  prior = prior,
  data = dat,
  sampler = sampler0,
  n_iter = 2e3,
  n_thin = 4,
  save = TRUE,
  prog = 0
) # samples for pseudo-target

(state2 <- mc_out$state)


### Additional prep

g_samples <- sapply(mc_out$sims, function(x) x$g)

if (isTRUE(logG)) {
  samples_use <- log(g_samples)
} else {
  samples_use <- g_samples
}

if (type %in% c("rw", "stepping", "latent")) {
  # will require tuning

  tune_bnds_init <- c(0.2 * sd(samples_use), 0.8 * diff(range(samples_use)))

  if (type == "latent") {
    state2$latent_s <- mean(tune_bnds_init)
    tune_bnds_init <- rev(1.0 / tune_bnds_init) # uses rate parameter
  }

  sampler_tuned$g$type <- type
  sampler_tuned$g$subtype <- subtype

  mc_tune <- tune(
    state = state2,
    prior = prior,
    data = dat,
    sampler = sampler_tuned,
    param = "g",
    bnds_init = tune_bnds_init,
    ess_log = ess_log,
    n_iter = 1e3,
    n_grid = 3,
    n_rep = 2,
    n_rounds = 5,
    range_frac = 0.67,
    verbose = TRUE
  )

  # plot(mc_tune$lvals, mc_tune$lesps)
  # abline(v = log(mc_tune$val_opt), lty = 2)

  sampler_tuned$g <- mc_tune$sampler$g

  ### burn-in again with new sampler (so all timing runs begin after iter 10,000)
  mc_out <- mcmc_gprior(
    state = state0,
    prior = prior,
    data = dat,
    sampler = sampler0, # let the slice sampler get to stationarity quickly
    n_iter = 2e3,
    save = FALSE,
    prog = 1000
  )

  (state1 <- mc_out$state)

  if (type == "latent") {
    state1$latent_s <- mean(tune_bnds_init)
  }

  mc_out <- mcmc_gprior(
    state = state1,
    prior = prior,
    data = dat,
    sampler = sampler_tuned, # let the tuned sampler do the rest
    n_iter = 8e3,
    save = FALSE,
    prog = 1000
  )

  (state2 <- mc_out$state)
} else if (grepl("samples", subtype)) {
  # tune with samples

  tmp_pseu <- pseudo_opt(
    samples = samples_use,
    type = "samples",
    family = "t",
    degf = c(1, 5),
    lb = ifelse(logG, -Inf, 0.0),
    ub = ifelse(logG, log(prior$g_max), prior$g_max),
    utility_type = "AUC",
    plot = FALSE
  )

  sampler_tuned$g$type <- type
  sampler_tuned$g$subtype <- subtype
  sampler_tuned$g$pseudo <- tmp_pseu$pseudo
  sampler_tuned$g$loc <- tmp_pseu$pseudo$params$loc
  sampler_tuned$g$sc <- tmp_pseu$pseudo$params$sc
  sampler_tuned$g$degf <- tmp_pseu$pseudo$params$degf
  sampler_tuned$g$txt <- tmp_pseu$pseudo$txt
} else if (grepl("Laplace_analytic", subtype)) {
  sampler_tuned$g$type <- type
  sampler_tuned$g$subtype <- subtype

  if (isTRUE(logG)) {
    sampler_tuned$g$degf <- 5

    if (grepl("wide", subtype)) {
      sampler_tuned$g$sc_adj <- 1.2
    } else {
      sampler_tuned$g$sc_adj <- 1.0
    }
  } else {
    sampler_tuned$g$degf <- 1

    if (grepl("wide", subtype)) {
      sampler_tuned$g$sc_adj <- 1.5
    } else {
      sampler_tuned$g$sc_adj <- 1.0
    }
  }
}

### timing run

sampler_tuned$g$doextra <- run_info$n_extra > 0 # should target evaluation include superfluous matrix computations?
if (isTRUE(sampler_tuned$g$doextra)) {
  sampler_tuned$g$n_extra <- run_info$n_extra
}

mc_time <- time_gprior(
  state = state2,
  prior = prior,
  data = dat,
  param = "g",
  sampler = sampler_tuned,
  n_iter = n_iter,
  ess_log = ess_log
)

mc_time$timing
(samp_p_sec <- mc_time$timing$EffSamp / mc_time$timing$userTime)

sampler_tuned$g$doextra <- FALSE


### diagnostics
draws_g <- mc_time$draws
(g_SE <- summary(as.mcmc(draws_g))$statistics["Time-series SE"])
(g_mn <- mean(draws_g))
(g_sd <- sd(draws_g))

draws_lg <- log(draws_g)
(lg_SE <- unname(summary(as.mcmc(draws_lg))$statistics["Time-series SE"]))
(lg_mn <- mean(draws_lg))
(lg_sd <- sd(draws_lg))

if (type %in% c("rw", "stepping", "latent")) {
  tune_param <- mc_tune$val_opt
} else {
  tune_param <- NA
}

if (type == "Qslice") {
  mc_out <- mcmc_gprior(
    state = state2,
    prior = prior,
    data = dat,
    sampler = sampler_tuned,
    n_iter = 2e3,
    save = TRUE,
    prog = 0
  ) # samples for pseudo-target

  draws_u <- sapply(mc_out$extras, function(x) x$u)
  (AUC <- auc(u = draws_u))
} else {
  AUC <- NA
}

tempDf <- data.frame(
  target = target,
  type = type,
  subtype = subtype,
  n_extra = run_info$n_extra,
  rep = run_info$rep,
  run_id = run_id,
  ii = ii,
  n_iter = n_iter,
  g_mn = g_mn,
  g_sd = g_sd,
  g_SE = g_SE,
  lg_mn = lg_mn,
  lg_sd = lg_sd,
  lg_SE = lg_SE,
  nEval = mc_time$timing$nEval,
  ESS = mc_time$timing$EffSamp,
  userTime = mc_time$timing$userTime,
  sysTime = mc_time$timing$sysTime,
  elapsedTime = mc_time$timing$elapsedTime,
  sampPsec = samp_p_sec,
  tuneParam = tune_param,
  auc = AUC
)

write.table(
  tempDf,
  file = paste0(
    "output/",
    "target",
    target,
    "_type",
    type,
    "_subtype",
    subtype,
    "_nextra",
    run_info$n_extra,
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

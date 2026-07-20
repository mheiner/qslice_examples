args <- commandArgs(trailingOnly = TRUE)
ii <- as.numeric(args[1]) # job id
dte <- as.numeric(args[2])

##### for testing
# ii <- 5205
# dte <- 260527
#####

source("0_setup.R")
source("functions/tune.R")

sessionInfo()

load(paste0("input/schedule_all_", dte, ".rda"))
head(sched, n = 14)

run_id <- job_order[ii]
(run_info <- sched[run_id, ])
target <- as.character(run_info$target)
type <- as.character(run_info$type)
subtype <- as.character(run_info$subtype)

source(paste0("0_setup_", target, ".R"))
trials <- trials[[target]]

n_iter <- 50e3
(ess_log <- grepl("gamma", target) & !grepl("log", target)) # always report ESS on log scale for gamma, inv-gamma (regardless of whether log(x) is sampled)


set.seed(dte + run_info$rep)
n_train <- 2e3
samples_train <- truth$q(runif(n_train))

if (type %in% c("rw", "latent", "stepping")) {
  tune_bnds_init <- c(0.2 * sd(samples_train), 1.2 * diff(range(samples_train)))

  set.seed(ii)
  (x_0 <- truth$q(runif(1)))
  state <- list(x = x_0)

  if (type == "latent") {
    s_init <- mean(tune_bnds_init)
    state$s <- s_init
    tune_bnds_init <- rev(1.0 / tune_bnds_init) # uses rate parameter
  }

  settings <- trials[[type]]
  settings$tune_param <- mean(tune_bnds_init)

  mc_tune <- tune(
    state = state,
    log_target = truth$ld,
    support = c(truth$lb, truth$ub),
    type = type,
    settings = settings,
    bnds_init = tune_bnds_init,
    ess_log = ess_log,
    n_iter = 1000,
    n_grid = 3,
    n_rep = 2,
    n_rounds = 5,
    range_frac = 0.8,
    verbose = TRUE
  )

  # plot(mc_tune$lvals, mc_tune$lesps)
  # abline(v = log(mc_tune$val_opt), lty = 2)

  settings$tune_param <- mc_tune$val_opt
  run_info$algo_descrip <- paste0("tuned: ", mc_tune$val_opt)
} else if (type == "gess") {
  settings <- trials[[type]][[subtype]][["pseudo"]]
} else {
  # pseudo-target

  if (
    isTRUE(
      run_info["subtype"] %in% c("MSW_samples", "AUC_samples", "MM_Cauchy")
    )
  ) {
    # samples-based

    source("1_pseudo_from_samples.R") # pseudo-target from samples
    run_info$algo_descrip <- trials[[type]][[subtype]]$algo_descrip
  }

  settings <- trials[[type]][[subtype]]$pseudo$pseu
}


set.seed(dte + ii)
(x_0 <- truth$q(runif(1)))
state <- list(x = x_0, s = ifelse(type == "latent", s_init, NA))

## pre run (warm-up for MKL libraries) / burn-in for latent slice
mcmc_time <- sampler_time_eval(
  type = type,
  state = state,
  n_iter = 1000,
  lf_func = truth$ld,
  support = c(truth$lb, truth$ub),
  settings = settings,
  ess_log = ess_log
)

(state <- mcmc_time$state)

## timing run
mcmc_time_out <- sampler_time_eval(
  type = type,
  state = state,
  n_iter = n_iter,
  lf_func = truth$ld,
  support = c(truth$lb, truth$ub),
  settings = settings,
  ess_log = ess_log
)

(thin <- min(10 * n_iter / mcmc_time_out$tbl$EffSamp, 100) |> round())
indx_thin <- seq(thin, n_iter, by = thin)
thinDraws <- mcmc_time_out$draws[indx_thin]


tempDf <- data.frame(
  target = target,
  type = type,
  subtype = subtype,
  rep = run_info$rep,
  run_id = run_id,
  ii = ii,
  algo_descrip = as.character(run_info$algo_descrip),
  n_iter = n_iter,
  x_0 = x_0,
  nEval = mcmc_time_out$tbl$nEval,
  ESS = mcmc_time_out$tbl$EffSamp,
  userTime = mcmc_time_out$tbl$userTime,
  sysTime = mcmc_time_out$tbl$sysTime,
  elapsedTime = mcmc_time_out$tbl$elapsedTime,
  sampPsec = mcmc_time_out$tbl$EffSamp / mcmc_time_out$tbl$userTime,
  ks_pval = ks.test(thinDraws, truth$p)$p.value
)

## save summary
write.table(
  tempDf,
  file = paste0(
    "output/target",
    target,
    "_type",
    type,
    "_subtype",
    subtype,
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

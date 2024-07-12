args <- commandArgs(trailingOnly = TRUE)

data_type <- args[1] # either "river_long" or "sim"
round <- args[2] # either "tune" or "test"
ii <- as.numeric(args[3]) # positive integer
dte <- as.numeric(args[4]) # positive integer

### or for manual interactive run
# data_type <- "sim"
# round <- "tune"
# ii <- 1
# dte <- 1
###

library("mcmcse") # multivariate effective sample size
library("dplyr") # for subsetting data
library("coda")
library("truncnorm") # truncated normal distribution
library("crch") # truncated t distribution
library("qslice")

sessionInfo()

### if jobs scheduled
load(paste0("data/schedule_data", data_type, "_", round, ".rda"))
(run_id <- job_order[ii])
(run_info <- sched[run_id,])

## or for manual interactive run
# run_id <- 1
# run_info <- data.frame(param = "normal", value = Inf, rep = 1, 
#                        sampler = "pits", sample_size = 50)
###

(sampler_type <- run_info$sampler)
(n <- run_info$sample_size)

(mesg <- paste(round, data_type, n, sampler_type, run_info$param, run_info$value, run_info$rep, ii, dte, sep = "_"))


# Import the relevant functions
source("scripts/ct_gibbs/dhr_gibbs.R")

## sample/burnin must be proportional across methods for times to be comparable (includes burn in)
if (round == "tune") {
  (burnin <- ifelse(grepl("latent|multivariate", sampler_type), 10e3, 2e3))
  (n_iter <- ifelse(grepl("latent|multivariate", sampler_type), 40e3, 8e3))
  n_chains <- 2
} else {
  (burnin <- ifelse(grepl("latent|multivariate", sampler_type), 10e3, 2e3))
  (n_iter <- ifelse(grepl("latent|multivariate", sampler_type), 100e3, 20e3))
  n_chains <- 2
}

p <- 3 # p = nfreq + 1
source("0_simulate.R")

set.seed(dte + ii)

if (sampler_type %in% c("Qslice", "pits", "independence_metropolis")) {
  sample_parameters$pseu_family <- run_info$param
  sample_parameters$pseu_params$degf <- run_info$value
} else if (run_info$param %in% names(sample_parameters)) {
  sample_parameters[run_info$param] <- run_info$value
}

mc_time <- system.time({
  fit <- dhr(parameters = sample_parameters, sampler = sampler_type, subtype = subtype,
             verbose = TRUE)
})
mc_time

source("scripts/convergence_analysis.R")
mc_diag <- convergence_analysis(fit, verbose = TRUE)

(esps <- mean(mc_diag$ess$theta) / mc_time["user.self"] |> unname())
(esps_all <- mean(c(mc_diag$ess$theta, mc_diag$ess$phase, mc_diag$ess$sig2)) / mc_time["user.self"] |> unname())

(mess <- mcmcse::multiESS(mc_diag$mc_list$theta_dg[[1]]))
(mesps <- mcmcse::multiESS(mc_diag$mc_list$theta_dg[[1]]) / mc_time["user.self"] |> unname())

mc_summary <- list()
mc_summary[["theta"]] <- summary(mc_diag$mc_list$theta_dg)
mc_summary[["phase"]] <- summary(mc_diag$mc_list$phases_dg)
mc_summary[["sig2"]] <- summary(mc_diag$mc_list$sigma2_dg)

rejections <- fit$rejections
gelman <- mc_diag$gelman
ess <- mc_diag$ess

save(mesg, round, run_info, run_id, ii, dte, 
     sample_parameters, data_sim, rejections, mc_time, gelman, 
     ess, esps, esps_all, mess, mesps,
     mc_summary,
     file = paste0("output/results_", mesg, ".rda"))

quit(save = "no")
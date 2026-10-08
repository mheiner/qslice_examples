# rm(list=ls())

#### inputs

## comment all targets out if running in a loop in 1_schedule_all.R
# target <- "normal"
# target <- "gamma"
# target <- "gammalog"
# target <- "igamma"
# target <- "igammalog"

n_rep <- 100 # number of replicate runs for each sampler/setting
wide_factor <- 4.0 # scale inflation for methods using a "diffuse" pseudo-target

#########

library("qslice")
library("tidyverse")

source(paste0("0_setup_", target, ".R"))

trials <- list()
competitors <- c("rw", "stepping", "latent")

for (cc in competitors) {
  trials[[cc]] <- list(target = target, type = cc, subtype = NA)
}

### create schedule
sched0 <- list()

init_lapprox <- 0.4 # initial value for optimization of Laplace approximation (can use across all targets)
source("1_setup_trials_pseudo.R") # methods with pseudo-targets won't require initial calibration

for (cc in competitors) {
  sched0[[cc]] <- expand.grid(
    target = target,
    type = cc,
    subtype = NA,
    stringsAsFactors = FALSE,
    algo_descrip = ""
  )
}

sched1 <- do.call(rbind, sched0)
rownames(sched1) <- NULL
sched1

sched <- lapply(1:n_rep, function(i) cbind(sched1, rep = i)) %>%
  do.call(rbind, .)
head(sched, n = 20)
tail(sched, n = 20)

(n_jobs <- nrow(sched))

save(
  file = paste0("input/schedule_target_", target, ".rda"),
  target,
  n_rep,
  trials,
  sched,
  n_jobs
)

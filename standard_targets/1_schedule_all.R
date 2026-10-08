rm(list = ls())

dte <- 261008
set.seed(dte)

# targets <- c("normal", "gamma", "igamma")
targets <- c("normal", "gamma", "igamma", "gammalog", "igammalog")

#### Optionally run here, in which case all target assignments should be commented out in 1_setup_trials.R
# for (target in targets) {
#   # can take time for pseudo-target optimization; optionally run the script separately on selected targets
#   source("1_setup_trials.R")
# }
# rm(list = setdiff(ls(), c("targets", "dte")))
####

sched_list <- list()
trials <- list()

for (tg in targets) {
  # requires 1_setup_trials.R be run for each target first
  tmp_env <- new.env()
  load(paste0("input/schedule_target_", tg, ".rda"), envir = tmp_env)
  sched_list[[tg]] <- tmp_env$sched
  trials[[tg]] <- tmp_env$trials
  rm(tmp_env)
}
rm(tg)

sched <- do.call(rbind, sched_list)

rm(list = setdiff(ls(), c("targets", "sched", "dte", "trials")))
ls()

str(sched)
rownames(sched) <- NULL
head(sched, n = 20)
tail(sched, n = 20)

(n_jobs <- nrow(sched))

job_order <- sample(n_jobs, size = n_jobs, replace = FALSE)

save(
  file = paste0("input/schedule_all_", dte, ".rda"),
  sched,
  n_jobs,
  job_order,
  trials
)

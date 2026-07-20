rm(list = ls())

dte <- 260704
set.seed(dte)

n_reps <- 50

Ks <- c(20)
targets <- c("unordered", "decreasing")
pseudo_types <- c("data", "samples")

samplers <- c(
     "QSS_stick",
     "QSS_logit",
     "IMH_stick",
     "IMH_logit",
     "GESS",
     "GPSS",
     "HRSS",
     "LSS",
     "RWM"
)

sched0 <- expand.grid(
     pseudo = pseudo_types,
     sampler = samplers,
     target = targets,
     K = Ks
)
dim(sched0)

sched <- do.call(rbind, lapply(1:n_reps, function(i) cbind(sched0, rep = i)))

str(sched)
head(sched, n = 40)
tail(sched, n = 40)

(n_jobs <- nrow(sched))

job_order <- sample(n_jobs, size = n_jobs, replace = FALSE)

save(
     file = paste0("input/schedule_all_", dte, ".rda"),
     sched,
     n_jobs,
     job_order
)

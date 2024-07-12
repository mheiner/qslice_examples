rm(list=ls())
library("tidyverse")

dte_sched <- 240603
dte_sched <- 2406033
dte_sched <- 240606
set.seed(dte_sched)

n_reps <- 20

# datatype <- "high-smooth"; sample_sizes <- 60
# datatype <- "low-smooth"
datatype <- "river"; sample_sizes <- 60
datatype <- "river_long"; sample_sizes <- 117
# datatype <- "river_long"; sample_sizes <- 190
# datatype <- "river_long"; sample_sizes <- 40
# sample_sizes <- c(60, 120, 240)

types <- "multivariate_slice"
types <- "latent_slice"
types <- c("Qslice", "pgibbs")
types <- c("pits", "independence_metropolis")
types <- c("ffbs", "independence_metropolis", "Qslice", "pgibbs")
types <- c("ffbs", "multivariate_slice", "latent_slice", "independence_metropolis", "Qslice", "pgibbs")
types <- c("multivariate_slice", "latent_slice", "independence_metropolis", "pits", "pgibbs")
types <- c("multivariate_slice", "latent_slice", "pits", "Qslice", "pgibbs")

tab_data <- expand.grid(data = datatype, sample_size = sample_sizes)
tab_data

tab_samplers_list <- list(
  multivariate_slice = data.frame(sampler = "multivariate_slice",
                                  param = "width",
                                  value = c(0.2, 0.5, 1.0, 2.0, 5.0)),
  latent_slice = data.frame(sampler = "latent_slice",
                            param = "latent_scale",
                            value = c(0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0)),
  independence_metropolis = data.frame(sampler = rep("independence_metropolis", 2),
                                       param = c(rep("t", 1), "normal"),
                                       value = c(5, Inf)),
  pits = data.frame(sampler = rep("pits", 2),
                      param = c(rep("t", 1), "normal"),
                      value = c(5, Inf)),
  Qslice = data.frame(sampler = rep("Qslice", 2),
                      param = c("t", "normal"),
                      value = c(5, Inf)),
  pgibbs = data.frame(sampler = "pgibbs",
                      param = "n_particles",
                      value = c(10, 15, 20, 25, 50))
)
tab_samplers_list
tab_samplers <- do.call(rbind, tab_samplers_list)
rownames(tab_samplers) <- NULL
tab_samplers

tab_experiments <- cross_join(tab_data, tab_samplers)
head(tab_experiments, n = 30); tail(tab_experiments, n = 30)

sched <- lapply(1:n_reps, function(i) cbind(tab_experiments, rep = i)) %>% do.call(rbind, .)
head(sched, n = 50); tail(sched, n = 50)

(n_jobs <- nrow(sched))

job_order <- sample(n_jobs, size = n_jobs, replace = FALSE)

save(file = paste0("data/schedule_data", datatype, "_tune.rda"),
     sched, n_jobs, job_order, dte_sched)

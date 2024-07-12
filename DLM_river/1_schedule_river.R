rm(list=ls())
library("tidyverse")

dte_sched <- 240604
set.seed(dte_sched)

n_reps <- 50

datatype <- "river_long"; sample_sizes <- 117
# datatype <- "sim"; sample_sizes <- 50

types <- c("ffbs", "multivariate_slice", "latent_slice", "independence_metropolis", "pits", "Qslice", "pgibbs")

tab_data <- expand.grid(data = datatype, sample_size = sample_sizes)
tab_data

tab_samplers_list <- list(
  ffbs = data.frame(sampler = "ffbs", param = "", value = 0),
  multivariate_slice = data.frame(sampler = "multivariate_slice",
                                  param = "width",
                                  value = c(1.5)),
  latent_slice = data.frame(sampler = "latent_slice",
                            param = "latent_scale",
                            value = c(0.3)),
  independence_metropolis = data.frame(sampler = rep("independence_metropolis", 2),
                                       param = c("t", "normal"),
                                       value = c(5, Inf)),
  pits = data.frame(sampler = rep("pits", 2),
                      param = c("t", "normal"),
                      value = c(5, Inf)),
  Qslice = data.frame(sampler = rep("Qslice", 2),
                      param = c("t", "normal"),
                      value = c(5, Inf)),
  pgibbs = data.frame(sampler = "pgibbs",
                      param = "n_particles",
                      value = c(25))
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

(outfile <- paste0("data/schedule_data", datatype, "_test.rda"))
save(file = outfile,
     sched, n_jobs, job_order, dte_sched)

rm(list=ls())

dte <- 240708

round <- "tune"
round <- "test"

data_types <- "river_long"

files_use <- list.files("output", pattern = paste0("results_", round, "_(", 
                                                   paste(data_types, collapse = "|"), 
                                                   ").*_", dte, ".rda"))
(nfiles <- length(files_use))

(files_sched <- list.files("data", pattern = paste0("schedule_data", 
                                                   paste(data_types, collapse = "|"), 
                                                   "_", round, ".rda")))

sched_list <- list()
for (i in 1:length(files_sched)) {
  load(paste0("data/", files_sched[i]))
  sched_list[[i]] <- sched
}

dat <- do.call(rbind, sched_list)
str(dat); head(dat)

(nrun <- nrow(dat))

dat$ii <- NA
dat$job <- NA
dat$n_iter_chain <- NA
dat$n_chains <- NA
dat$sec <- NA
dat$ess_theta_mean <- NA
dat$ess_theta_min <- NA
dat$mess <- NA
dat$rejections_mean <- NA
dat$GelmanU_theta_mean <- NA
dat$GelmanU_theta_max <- NA
dat$theta1post_mean <- NA
dat$theta2post_mean <- NA
dat$theta3post_mean <- NA
dat$theta4post_mean <- NA
dat$theta5post_mean <- NA
dat$theta6post_mean <- NA
dat$theta7post_mean <- NA
dat$theta8post_mean <- NA
dat$theta9post_mean <- NA

dat$theta1post_se <- NA
dat$theta2post_se <- NA
dat$theta3post_se <- NA

dat$theta10post_mean <- NA
dat$theta10post_se <- NA
dat$theta11post_mean <- NA
dat$theta11post_se <- NA
dat$theta12post_mean <- NA
dat$theta12post_se <- NA
dat$theta13post_mean <- NA
dat$theta13post_se <- NA
dat$theta14post_mean <- NA
dat$theta14post_se <- NA
dat$theta15post_mean <- NA
dat$theta15post_se <- NA
dat$theta16post_mean <- NA
dat$theta16post_se <- NA

for (i in 1:nfiles) {

  load(paste0("output/", files_use[i]))  
  
  indx_now <- run_id
  
  if (length(indx_now) != 1) {
    cat("WARNING: not exactly one row match.\n")
    cat("indx:", indx_now, "\n")
    cat("iter:", i, "\n")
    print(run_info)
  }

  dat[indx_now, "ii"] <- ii
  dat[indx_now, "job"] <- job_order[ii]
  dat[indx_now, "n_iter_chain"] <- sample_parameters$n_iter
  dat[indx_now, "n_chains"] <- sample_parameters$n_chains
  dat[indx_now, "sec"] <- mc_time["user.self"] |> unname()
  dat[indx_now, "ess_theta_mean"] <- mean(ess$theta)
  dat[indx_now, "ess_theta_min"] <- min(ess$theta)
  dat[indx_now, "mess"] <- mess
  dat[indx_now, "rejections_mean"] <- mean(rejections)
  dat[indx_now, "GelmanU_theta_mean"] <- mean(gelman$theta$psrf[,"Upper C.I."])
  dat[indx_now, "GelmanU_theta_max"] <- max(gelman$theta$psrf[,"Upper C.I."])
  dat[indx_now, "theta1post_mean"] <- mc_summary$theta$statistics[1, "Mean"]
  dat[indx_now, "theta2post_mean"] <- mc_summary$theta$statistics[2, "Mean"]
  dat[indx_now, "theta3post_mean"] <- mc_summary$theta$statistics[3, "Mean"]
  dat[indx_now, "theta4post_mean"] <- mc_summary$theta$statistics[4, "Mean"]
  dat[indx_now, "theta5post_mean"] <- mc_summary$theta$statistics[5, "Mean"]
  dat[indx_now, "theta6post_mean"] <- mc_summary$theta$statistics[6, "Mean"]
  dat[indx_now, "theta7post_mean"] <- mc_summary$theta$statistics[7, "Mean"]
  dat[indx_now, "theta8post_mean"] <- mc_summary$theta$statistics[8, "Mean"]
  dat[indx_now, "theta9post_mean"] <- mc_summary$theta$statistics[9, "Mean"]
  dat[indx_now, "theta1post_se"] <- mc_summary$theta$statistics[1, "Time-series SE"]
  dat[indx_now, "theta2post_se"] <- mc_summary$theta$statistics[2, "Time-series SE"]
  dat[indx_now, "theta3post_se"] <- mc_summary$theta$statistics[3, "Time-series SE"]
  dat[indx_now, "theta10post_mean"] <- mc_summary$theta$statistics[10, "Mean"]
  dat[indx_now, "theta10post_se"] <- mc_summary$theta$statistics[10, "Time-series SE"]
  dat[indx_now, "theta11post_mean"] <- mc_summary$theta$statistics[15, "Mean"]
  dat[indx_now, "theta11post_se"] <- mc_summary$theta$statistics[15, "Time-series SE"]
  dat[indx_now, "theta12post_mean"] <- mc_summary$theta$statistics[30, "Mean"]
  dat[indx_now, "theta12post_se"] <- mc_summary$theta$statistics[30, "Time-series SE"]
  dat[indx_now, "theta13post_mean"] <- mc_summary$theta$statistics[40, "Mean"]
  dat[indx_now, "theta13post_se"] <- mc_summary$theta$statistics[40, "Time-series SE"]
  dat[indx_now, "theta14post_mean"] <- mc_summary$theta$statistics[100, "Mean"]
  dat[indx_now, "theta14post_se"] <- mc_summary$theta$statistics[100, "Time-series SE"]
  dat[indx_now, "theta15post_mean"] <- mc_summary$theta$statistics[111, "Mean"]
  dat[indx_now, "theta15post_se"] <- mc_summary$theta$statistics[111, "Time-series SE"]
  dat[indx_now, "theta16post_mean"] <- mc_summary$theta$statistics[141, "Mean"]
  dat[indx_now, "theta16post_se"] <- mc_summary$theta$statistics[141, "Time-series SE"]
  
  cat(i, "of", nfiles, "\r")
}

head(dat); tail(dat)

dat$n_iter_tot <- dat$n_iter_chain * dat$n_chains

(outfile <- paste0("output/summary_", round, "_", paste(data_types, collapse = "-"), "_", dte, ".rda"))
save(dat, file = outfile)

(files_summary <- list.files("output", pattern = "summary"))


rm(list=ls())

(files_summary <- list.files("output", pattern = "summary"))



## interactive checking

sampler_now <- "ffbs"
sampler_now <- "multivariate_slice"
sampler_now <- "latent_slice"
sampler_now <- "independence_metropolis"
sampler_now <- "Qslice"
sampler_now <- "pgibbs"

Rhat_cutoff <- 1.05
Rhat_cutoff <- 1.5

ggplot(dat %>% filter(sampler == sampler_now),
       aes(x = n_iter_tot / sec, 
           y = paste(data, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = sampler_now) + ylab("")

ggplot(dat %>% filter(sampler == sampler_now),
       aes(x = GelmanU_theta_mean, 
           y = paste(data, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = sampler_now) + ylab("") + geom_vline(xintercept = Rhat_cutoff) +
  xlim(c(0.95, 3.0))

ggplot(dat %>% filter(sampler == sampler_now),
       aes(x = GelmanU_theta_max, 
           y = paste(data, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = sampler_now) + ylab("") + geom_vline(xintercept = Rhat_cutoff) +
  xlim(c(0.95, 3.0))

ggplot(dat %>% filter(sampler == sampler_now, GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_mean / sec, 
           y = paste(data, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = sampler_now) + ylab("")

ggplot(dat %>% filter(sampler == sampler_now, GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_min / sec, 
           y = paste(data, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = sampler_now) + ylab("")

ggplot(dat %>% filter(sampler == sampler_now, GelmanU_theta_mean < Rhat_cutoff),
       aes(x = mess * n_chains / sec, # mess is reported for only one chain 
           y = paste(data, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = sampler_now) + ylab("")

ggplot(dat %>% filter(sampler == sampler_now, GelmanU_theta_mean < Rhat_cutoff),
       aes(x = rejections_mean / n_iter_tot, 
           y = paste(data, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = sampler_now) + ylab("")


ggplot(dat %>% filter(data == "low-smooth" & sampler != "ffbs", GelmanU_theta_mean < Rhat_cutoff, sample_size == 60),
       aes(x = theta1post_mean, 
           y = paste(data, sampler, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = "estimation") + ylab("")

ggplot(dat %>% filter(data == "high-smooth", GelmanU_theta_mean < Rhat_cutoff, sample_size == 60),
       aes(x = theta4post_mean, 
           y = paste(data, sampler, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = "estimation") + ylab("")

ggplot(dat %>% filter(data == "river", GelmanU_theta_mean < Rhat_cutoff, sample_size == 60),
       aes(x = theta4post_mean, 
           y = paste(data, sampler, sprintf("%03d", sample_size), param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = "estimation") + ylab("")

dat %>% filter(data == "high-smooth", GelmanU_theta_mean < Rhat_cutoff, sample_size == 60) %>%
  group_by(data, sampler, param, value) %>% summarize(mean_th = mean(theta4post_mean), mean_se = mean(theta4post_se))


ggplot(dat %>% filter(value %in% c(0, 0.2, 5, Inf, 25), GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_mean / sec, 
           y = paste(data, sprintf("%03d", sample_size), sampler, param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = "Comparison") + ylab("")

ggplot(dat %>% filter(value %in% c(0, 0.2, 5, Inf, 25), GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_min / sec, 
           y = paste(data, sprintf("%03d", sample_size), sampler, param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = "Comparison") + ylab("")

ggplot(dat %>% filter(value %in% c(0, 0.2, 5, Inf, 25), GelmanU_theta_mean < Rhat_cutoff),
       aes(x = mess * n_chains / sec, 
           y = paste(data, sprintf("%03d", sample_size), sampler, param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = "Comparison") + ylab("")

ggplot(dat %>% filter(value %in% c(0, 0.2, 5, Inf, 25), GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_mean / (n_iter_chain * n_chains), 
           y = paste(data, sprintf("%03d", sample_size), sampler, param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = "Comparison") + ylab("")

ggplot(dat %>% filter(value %in% c(0, 0.2, 5, Inf, 25), GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_min / (n_iter_chain * n_chains), 
           y = paste(data, sprintf("%03d", sample_size), sampler, param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = "Comparison") + ylab("")

ggplot(dat %>% filter(value %in% c(0, 0.2, 5, Inf, 25), GelmanU_theta_mean < Rhat_cutoff),
       aes(x = mess / (n_iter_chain), 
           y = paste(data, sprintf("%03d", sample_size), sampler, param, sprintf("%03.02f", value)), 
           color = data)) +
  geom_point() + ggtitle(label = "Comparison") + ylab("")

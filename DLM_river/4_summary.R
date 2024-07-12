rm(list=ls())

library("tidyverse")

dte <- 240708 # R version 4.4.1

# round <- "tune"
round <- "test"

data_types <- "river_long"

(filename <- paste0("output/summary_", round, "_", paste(data_types, collapse = "-"), "_", dte, ".rda"))
load(filename)

dim(dat)
head(dat)

old <- options(pillar.sigfig = 5)
dat %>% group_by(sampler, param, value) %>%
  summarize(
    iter_p_sec_mn = mean(n_chains*1.1*n_iter_chain/sec), # includes burn-in that was 10% as many as iters kept
    iter_p_sec_sd = sd(n_chains*1.1*n_iter_chain/sec),
    gelmanUtheta_mean_mn = mean(GelmanU_theta_mean), 
    gelmanUtheta_mean_sd = sd(GelmanU_theta_mean),
    esps_mean_mn = mean(ess_theta_mean / sec),
    esps_mean_sd = sd(ess_theta_mean / sec),
    esps_min_mn = mean(ess_theta_min / sec),
    esps_min_sd = sd(ess_theta_min / sec),
    neval_p_iter_mn = mean((rejections_mean*n_chains + 2*2*1.1*n_iter_chain*n_chains) / (1.1*n_iter_chain*n_chains*2)), # add back 2 baseline evals assuming +10% burnin; the final *2 is because there are two parameter vectors
    neval_p_iter_sd = sd((rejections_mean*n_chains + 2*2*1.1*n_iter_chain*n_chains) / (1.1*n_iter_chain*n_chains*2))
  )
## NOTE: rejections_mean is averaging total rejections for both alpha vectors (burn & regular samples) over the two chains




ls()
head(dat); tail(dat)

Rhat_cutoff <- 1.005
Rhat_cutoff <- Inf

ggplot(dat,
       aes(x = GelmanU_theta_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Rhat mean") + ylab("") + 
  geom_vline(xintercept = Rhat_cutoff)

ggplot(dat,
       aes(x = GelmanU_theta_max, 
           y = paste(sampler, value),
           color = sampler)) +
  geom_point() + ggtitle(label = "Rhat max") + ylab("") + 
  geom_vline(xintercept = Rhat_cutoff)

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_mean / sec, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "ESS / sec") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_min / sec, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "ESS min / sec") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = mess * n_chains / sec, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "MESS / sec") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_mean / (n_iter_chain * n_chains), 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Effeciency") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = ess_theta_min / (n_iter_chain * n_chains), 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Efficiency (min)") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = mess / (n_iter_chain), 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Effeciency (multivariate)") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = rejections_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Rejections") + ylab("")




Rhat_cutoff <- 1.001
Rhat_cutoff <- 1.05

dat <- dat %>% filter(sampler != "ffbs")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta1post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 1") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta2post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 2") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta3post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 3") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta4post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 4") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta5post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 5") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta6post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 6") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta7post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 7") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta8post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 8") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta9post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 9") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta10post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 10") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta11post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 15") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta12post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 30") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta13post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 40") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta14post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 100") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta15post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 111") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta16post_mean, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation: theta 141") + ylab("")





ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta1post_se, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation SE: theta 1") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta2post_se, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation SE: theta 2") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta3post_se, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation SE: theta 3") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta10post_se, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation SE: theta 10") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta13post_se, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation SE: theta 40") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta14post_se, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation SE: theta 100") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta15post_se, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation SE: theta 111") + ylab("")

ggplot(dat %>% filter(GelmanU_theta_mean < Rhat_cutoff),
       aes(x = theta16post_se, 
           y = paste(sampler, value), 
           color = sampler)) +
  geom_point() + ggtitle(label = "Estimation SE: theta 141") + ylab("")


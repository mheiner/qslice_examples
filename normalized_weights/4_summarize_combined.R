rm(list = ls())
library("tidyverse")

dte <- 260704 # f06; K = 20; C ~ lognormal(0, sig = 2.0); burn-in with GESS

# load(paste0("input/schedule_all_", dte, ".rda"))
dat0 <- read.csv(paste0("output/combined_all_", dte, ".csv"))
dat0$sampler <- factor(
  dat0$sampler,
  levels = c(
    "RWM",
    "IMH_stick",
    "IMH_logit",
    "HRSS",
    "LSS",
    "GESS",
    "GPSS",
    "QSS_stick",
    "QSS_logit"
  )
)
dat0$target <- factor(dat0$target, levels = c("unordered", "decreasing"))
dat0 <- dat0 %>% arrange(run_id)

head(dat0)
str(dat0)

table(dat0$sampler, dat0$target, dat0$pseudo, dat0$K)


## select
targets_now <- "unordered"
targets_now <- "decreasing"
targets_now <- c("unordered", "decreasing")
pseudos_now <- c("data", "samples")
pseudos_now <- "data"
pseudos_now <- "samples"
K_now <- 10
K_now <- 20
K_now <- 25

dat <- dat0 %>%
  filter(target %in% targets_now, K %in% K_now, pseudo %in% pseudos_now)
dat$pseudo_alpha <- ifelse(dat$pseudo == "samples", 0.9, 1.0)

dim(dat)


options(pillar.sigfig = 5)


## summarize tuning parameters
ggplot(
  dat %>% filter(sampler %in% c("RWM", "HRSS", "LSS", "GPSS")),
  aes(x = tune_param_val, y = sampler, color = pseudo)
) +
  geom_jitter(height = 0.1, width = 0, size = 0.2) +
  # facet_grid(~ type, scale = "free_x") +
  theme_bw()


## number of evauations
plt <- ggplot(
  dat,
  aes(
    x = nEval_mean,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  xlab("Evaluations per iter") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt

dat %>%
  group_by(sampler) %>%
  summarize(neval_mn = mean(nEval_mean), neval_sd = sd(nEval_mean))


## effective samples per second
plt <- ggplot(
  dat,
  aes(
    x = ESSZ_mean_total / userTime_total,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("Effective samples per second") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt


## min effective samples per second
plt <- ggplot(
  dat,
  aes(
    x = ESpSZ_min_avg,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("Min effective samples per second") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt


## R hat
plt <- ggplot(
  dat %>% filter(Rhat_u95_avg < 1.05),
  aes(
    x = Rhat_u95_avg,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("R hat") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt

dat %>%
  group_by(sampler) %>%
  summarize(
    Rhat_mn = mean((Rhat_u95_avg - 1) * 1000, na.rm = TRUE),
    Rhat_sd = sd((Rhat_u95_avg - 1) * 1000, na.rm = TRUE)
  )
table(dat$target, dat$sampler, is.na(dat$Rhat_u95_avg))


## effective samples per sample
plt <- ggplot(
  dat,
  aes(
    x = espitZ_avg,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("Effective samples per sample") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt

## IAT
plt <- ggplot(
  dat,
  aes(
    x = IATZ_avg,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("IAT") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt

## RMSE
plt <- ggplot(
  dat,
  aes(
    x = RMSE_logits_avg,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("RMSE") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt


### now with respect to reference

ref_sampler <- "GESS"

dat_ref <- dat %>% filter(sampler %in% ref_sampler)
dat <- left_join(
  dat,
  dat_ref,
  by = c("target", "K", "pseudo", "rep"),
  suffix = c("", "_ref")
)

head(dat, n = 50)


## number of evauations
plt <- ggplot(
  dat,
  aes(
    x = nEval_mean - nEval_mean_ref,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  xlab("Evaluations per iter (difference)") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt


## effective samples per second
plt <- ggplot(
  dat,
  aes(
    x = (ESSZ_mean_total / userTime_total) /
      (ESSZ_mean_total_ref / userTime_total_ref),
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("Effective samples per second (ratio)") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt

# ggsave(filename = paste0("plots/esps_rat_", dte, ".pdf"), width = 7, height = 8)

dat %>%
  group_by(sampler) %>%
  summarize(
    ESpSr_mn = mean(
      (ESSZ_mean_total / userTime_total) /
        (ESSZ_mean_total_ref / userTime_total_ref)
    ),
    ESpSr_sd = sd(
      (ESSZ_mean_total / userTime_total) /
        (ESSZ_mean_total_ref / userTime_total_ref)
    )
  )


## min effective samples per second
plt <- ggplot(
  dat,
  aes(
    x = ESpSZ_min_avg / ESpSZ_min_avg_ref,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("Min effective samples per second (ratio)") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt

dat %>%
  group_by(sampler) %>%
  summarize(
    ESpSminR_mn = mean(ESpSZ_min_avg / ESpSZ_min_avg_ref),
    ESpSminR_sd = sd(ESpSZ_min_avg / ESpSZ_min_avg_ref)
  )


## R hat
plt <- ggplot(
  dat %>% filter(Rhat_u95_avg < 1.1),
  aes(
    x = Rhat_u95_avg / Rhat_u95_avg_ref,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("R hat (ratio)") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt

## effective samples per sample
plt <- ggplot(
  dat,
  aes(
    x = espitZ_avg / espitZ_avg_ref,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("Effective samples per sample (ratio)") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt

dat %>%
  group_by(sampler) %>%
  summarize(
    ESpIt_mn = mean(espitZ_avg / espitZ_avg_ref),
    ESpIt_sd = sd(espitZ_avg / espitZ_avg_ref)
  )


## RMSE
plt <- ggplot(
  dat %>% filter(RMSE_logits_avg / RMSE_logits_avg_ref < 1.5),
  aes(
    x = RMSE_logits_avg / RMSE_logits_avg_ref,
    y = sampler,
    fill = sampler,
    color = pseudo,
    alpha = pseudo_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + # theme(legend.position = "none") +
  scale_color_manual(values = c("gray20", "gray"), drop = FALSE) +
  scale_alpha(guide = "none") +
  xlab("RMSE (ratio)") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target, scales = "free_x")

plt

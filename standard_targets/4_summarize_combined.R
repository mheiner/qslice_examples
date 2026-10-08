rm(list = ls())
library("tidyverse")

targets <- "all"
# targets <- c("normal", "gamma", "igamma")
# targets <- c("gamma", "gammalog", "igamma", "igammalog")

dte <- 260527 # f06 server, 32 jobs in parallel
dte <- 260813 # f06 server, 32 jobs in parallel
dte <- 260822 # f06 server, 24 jobs in parallel, error handling
dte <- 260928 # f07 server, 32 jobs in parallel, diagnostic mode off
dte <- 261006 # f06 server, 32 jobs parallel, diagnostic mode off
dte <- 261008 # f04 server, 32 jobs parallel, diagnostic mode off; multiple psuedo_opt tries


if (length(targets) == 1 && targets == "all") {
  dat <- read.csv(paste0("output/combined_all_", dte, ".csv"))
} else {
  datl <- list()
  for (tg in targets) {
    datl[[tg]] <- read.csv(paste0(
      "output/combined_target",
      tg,
      "_",
      dte,
      ".csv"
    ))
  }
  dat <- do.call(rbind, datl)
}

str(dat)
hist(dat$ks_pval, breaks = 30)

## summarize Kolmogorov-Smirnov test statistics
dat$type_sub <- paste(dat$type, dat$subtype)
ks_rej <- dat %>%
  group_by(type_sub, target) %>%
  summarize(ks_rej = mean(ks_pval < 0.05), n = n())
ggplot(ks_rej, aes(x = jitter(ks_rej), y = type_sub, col = target)) +
  geom_point()

summary(ks_rej)


## summarize tuning parameters
ggplot(
  dat %>%
    filter(
      target %in% c("normal", "gamma", "igamma"),
      type %in% c("rw", "stepping", "latent")
    ) %>%
    mutate(tune = as.numeric(gsub("tuned: ", "", algo_descrip))),
  aes(x = tune, y = target)
) +
  geom_jitter(height = 0.1, width = 0, size = 0.2) +
  facet_grid(~type, scale = "free_x") +
  theme_bw()


## summarize algorithm settings
algo_descrips <- dat %>%
  filter(!grepl("samples|MM", subtype)) %>%
  group_by(type_sub, target) %>%
  summarize(n = n(), t = length(unique(algo_descrip)))
print(algo_descrips, n = 50)
dat %>%
  filter(
    !grepl("samples|MM", subtype),
    target %in% c("normal", "gamma", "igamma")
  ) %>%
  group_by(type_sub, target) %>%
  summarize(n = n(), t = length(unique(algo_descrip))) %>%
  print(., n = 50)

dat <- dat %>%
  mutate(algo = paste(type, subtype), algo_all = paste(type, algo_descrip))
dat$algo <- gsub(" NA", "", x = dat$algo)
unique(dat$algo)
unique(dat$algo_al)


dat$algoF <- factor(
  dat$algo,
  levels = rev(c(
    "rw",
    "imh AUC",
    "imh AUC_wide",
    "stepping",
    "gess AUC",
    "latent",
    "Qslice MSW",
    "Qslice AUC",
    "Qslice MSW_samples",
    "Qslice AUC_samples",
    "Qslice AUC_wide",
    "Qslice Laplace_Cauchy",
    "Qslice MM_Cauchy"
  )),
  labels = rev(c(
    "Random walk",
    "Independence M-H: AUC",
    "Independence M-H: AUC-diffuse",
    "Stepping out & shrinkage (Neal, 2003)",
    "Generalized elliptical (Nishihara et al., 2014)",
    "Latent slice (Li and Walker, 2023)",
    "Quantile slice: MSW",
    "Quantile Slice: AUC",
    "Quantile slice: MSW-samples",
    "Quantile Slice: AUC-samples",
    "Quantile slice: AUC-diffuse",
    "Quantile slice: Laplace-Cauchy",
    "Quantile Slice: MM-Cauchy"
  ))
)


dat$typeF <- case_match(
  dat$type,
  c("gess", "latent", "stepping") ~ "Slice",
  "imh" ~ "IMH",
  "Qslice" ~ "Quantile slice",
  "rw" ~ "Rand walk"
) %>%
  factor(., levels = c("Rand walk", "Slice", "Quantile slice", "IMH"))

dat$targetlab <- case_match(
  dat$target,
  "normal" ~ "Normal target",
  "gamma" ~ "Gamma target",
  "gammalog" ~ "Gamma-log",
  "igamma" ~ "Inverse-gamma target",
  "igammalog" ~ "Inverse Gamma-log"
)

dat$target_base <- case_match(
  dat$target,
  "normal" ~ "Normal target",
  c("gamma", "gammalog") ~ "Gamma target",
  c("igamma", "igammalog") ~ "Inverse-gamma target"
) %>%
  factor(., levels = c("Normal target", "Gamma target", "Inverse-gamma target"))

dat$target_tx <- ifelse(
  grepl("log", dat$target),
  "Log transform",
  "Original"
) %>%
  as.factor()
dat$target_tx_alpha <- ifelse(dat$target_tx == "Log transform", 0.5, 1.0)


plt <- ggplot(
  dat %>% filter(target %in% c("normal", "gamma", "igamma")),
  aes(x = sampPsec / 1e3, y = algoF, fill = typeF),
  color = "gray"
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  xlab("Effective samples (in thousands) per second") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target_base, scales = "free_x")

plt

ggsave(
  plot = plt,
  width = 9,
  height = 5.5,
  filename = paste0("plots/ESPS_primary_", dte, ".pdf")
)


plt <- ggplot(
  dat %>% filter(target %in% c("normal", "gamma", "igamma")),
  aes(x = ESS / n_iter, y = algoF, fill = typeF),
  color = "gray"
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  xlab("Effective samples (in thousands) per iteration") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target_base, scales = "free_x")

plt

ggsave(
  plot = plt,
  width = 9,
  height = 5.5,
  filename = paste0("plots/ESpIt_primary_", dte, ".pdf")
)


plt <- ggplot(
  dat %>% filter(target %in% c("normal", "gamma", "igamma")),
  aes(x = sampPsec / 1e3, y = algoF, color = ii)
) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  xlab("Effective samples (in thousands) per second") +
  ylab("") +
  facet_wrap(~target_base, scales = "free_x")

plt


## investigate issue with Qslice AUC-samples on inverse-gamma
dat_ig_auc <- dat %>%
  filter(target == "igamma", algo == "Qslice AUC_samples") %>%
  select(ii, target, type, subtype, algo_descrip, nEval, sampPsec) %>%
  arrange(sampPsec)
dat_ig_auc$loc <- as.numeric(sub(
  ".*loc = ([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?).*",
  "\\1",
  dat_ig_auc$algo_descrip
))
dat_ig_auc$sc <- as.numeric(sub(
  ".*sc = ([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?).*",
  "\\1",
  dat_ig_auc$algo_descrip
))
ggplot(dat_ig_auc, aes(x = sampPsec, y = nEval)) + geom_point()
ggplot(dat_ig_auc, aes(x = sampPsec, y = loc)) + geom_point()
ggplot(dat_ig_auc, aes(x = sampPsec, y = loc, col = sc)) + geom_point()

## investigate issue with Qslice MSW on gamma
dat_ga_msw <- dat %>%
  filter(target %in% c("gamma", "gammalog"), grepl("Qslice MSW", algo)) %>%
  select(ii, target, type, subtype, algo, algo_descrip, nEval, sampPsec) %>%
  arrange(sampPsec)
dat_ga_msw$loc <- as.numeric(sub(
  ".*loc = ([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?).*",
  "\\1",
  dat_ga_msw$algo_descrip
))
dat_ga_msw$sc <- as.numeric(sub(
  ".*sc = ([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?).*",
  "\\1",
  dat_ga_msw$algo_descrip
))
ggplot(
  dat_ga_msw %>% filter(grepl("log", target), grepl("samples", algo)),
  aes(x = loc, y = sc)
) +
  geom_point()


plt <- ggplot(
  dat %>% filter(target %in% c("gamma", "igamma", "gammalog", "igammalog")),
  aes(
    x = sampPsec / 1e3,
    y = algoF,
    fill = typeF,
    color = target_tx,
    alpha = target_tx_alpha
  )
) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  xlab("Effective samples (in thousands) per second") +
  ylab("") +
  labs(color = "", fill = "") +
  facet_wrap(~target_base, scales = "fixed") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(
  plot = plt,
  width = 9,
  height = 8,
  filename = paste0("plots/ESPS_gammas_w_tnx_", dte, ".pdf")
)


## compare performance with MSW and AUC

dat_q <- dat %>% filter(type == "Qslice")

str(dat_q)
head(dat_q)

dat_q$fam <- substr(dat_q$algo_descrip, 1, 1)
head(dat_q$fam, n = 200)

dat_q$loc <- as.numeric(sub(
  ".*loc = ([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?).*",
  "\\1",
  dat_q$algo_descrip
))
dat_q$sc <- as.numeric(sub(
  ".*sc = ([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?).*",
  "\\1",
  dat_q$algo_descrip
))
dat_q$df <- NA
dat_q$df[which(dat_q$fam == "C")] <- 1
dat_q$df[which(dat_q$fam == "t")] <- as.numeric(sub(
  ".*df = ([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?).*",
  "\\1",
  dat_q$algo_descrip[which(dat_q$fam == "t")]
))
table(dat_q$df)

dat_q$lb <- ifelse(dat_q$target %in% c("gamma", "igamma"), 0, -Inf)

head(dat_q)

dat_q_unique <- dat_q %>%
  select(target, type, subtype, algo_descrip, fam, loc, sc, df, lb) %>%
  distinct()

str(dat_q_unique)
head(dat_q_unique)

table(dat_q_unique$target, dat_q_unique$subtype)

library("qslice")

dat_q_unique$AUC <- NA
dat_q_unique$MSW <- NA

pb <- utils::txtProgressBar(min = 0, max = nrow(dat_q_unique), style = 3)

for (i in 1:nrow(dat_q_unique)) {
  target_now <- dat_q_unique[i, "target"]
  source(paste0("0_setup_", target_now, ".R"))
  pseudo_now <- qslice::pseudo_list(
    family = "t",
    params = list(
      loc = dat_q_unique[i, "loc"],
      sc = dat_q_unique[i, "sc"],
      df = dat_q_unique[i, "df"]
    ),
    lb = dat_q_unique[i, "lb"],
    ub = Inf
  )
  dat_q_unique[i, "MSW"] <- qslice::utility_pseudo(
    pseudo = pseudo_now,
    log_target = truth$ld,
    type = "function",
    utility_type = "MSW",
    plot = FALSE
  )
  dat_q_unique[i, "AUC"] <- qslice::utility_pseudo(
    pseudo = pseudo_now,
    log_target = truth$ld,
    type = "function",
    utility_type = "AUC",
    plot = FALSE
  )

  utils::setTxtProgressBar(pb, i)
}

close(pb)

head(dat_q_unique)
tail(dat_q_unique)

dat_q_metric <- left_join(
  dat_q,
  dat_q_unique,
  by = join_by(target, type, subtype, algo_descrip)
)

head(dat_q_metric)
tail(dat_q_metric)

dat_q_metric_show <- dat_q_metric %>%
  filter(target %in% c("normal", "gamma", "igamma"))

dat_q_metric_show$samples <- grepl("samples|MM", dat_q_metric_show$subtype)

p_auc <- ggplot(
  dat_q_metric_show,
  aes(x = AUC, y = sampPsec / 1e3, col = target, pch = samples)
) +
  geom_jitter(width = 0.01, height = 0, size = 0.6) +
  theme_bw() +
  ylab("ES/sec (thousands)") +
  xlim(-0.01, 1.01)

# ggsave(filename = paste0("plots/ESpSvAUC_", dte, ".pdf"), width = 5, height = 4)

p_msw <- ggplot(
  dat_q_metric_show,
  aes(x = MSW, y = sampPsec / 1e3, col = target, pch = samples)
) +
  geom_jitter(width = 0.01, height = 0, size = 0.6) +
  theme_bw() +
  ylab("ES/sec (thousands)") +
  xlim(-0.01, 1.01) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3))
  )

library("cowplot")
leg <- get_legend(p_msw + theme(legend.box.margin = margin(0, 0, 0, 0)))

cowplot::plot_grid(
  p_auc + theme(legend.position = "none"),
  p_msw + theme(legend.position = "none"),
  leg,
  rel_widths = c(1, 1, 0.3),
  nrow = 1
)

# ggsave(
#   filename = paste0("plots/ESpSvAUCMSW_", dte, ".pdf"),
#   width = 8,
#   height = 4
# )

ggplot(
  dat_q_metric_show,
  aes(x = AUC, y = ESS / n_iter, col = target, pch = samples)
) +
  geom_jitter(width = 0.01, height = 0) +
  theme_bw() +
  xlim(0, 1)

ggplot(
  dat_q_metric_show,
  aes(x = MSW, y = ESS / n_iter, col = target, pch = samples)
) +
  geom_jitter(width = 0.01, height = 0) +
  theme_bw() +
  xlim(0, 1)

p_aucmsw <- ggplot(
  dat_q_metric_show,
  aes(x = MSW, y = AUC, col = target, pch = samples)
) +
  geom_jitter(width = 0.01, height = 0.01, size = 0.6) +
  theme_bw() +
  xlim(-0.01, 1.01) +
  ylim(-0.01, 1.01) +
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed")


cowplot::plot_grid(
  p_msw + theme(legend.position = "none"),
  p_auc + theme(legend.position = "none"),
  p_aucmsw + theme(legend.position = "none"),
  leg,
  rel_widths = c(1, 1, 1, 0.5),
  nrow = 2,
  ncol = 2
)

ggsave(
  filename = paste0("plots/ESpSvAUCMSW3_", dte, ".pdf"),
  width = 8,
  height = 6
)

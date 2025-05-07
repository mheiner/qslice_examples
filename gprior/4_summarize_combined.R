rm(list=ls()); dev.off()
library("tidyverse")

targets <- "all"

dte <- 250502


if (targets == "all") {
  dat <- read.csv(paste0("output/combined_all_", dte, ".csv"))
} else {
  datl <- list()
  for (tg in targets) {
    datl[[tg]] <- read.csv(paste0("output/combined_target", tg, "_", dte, ".csv"))
  }
  dat <- do.call(rbind, datl)
}

str(dat)
setdiff(1:max(dat$run_id), dat$run_id)


## summarize algorithm settings
dat <- dat %>% mutate(algo = paste(type, subtype))
dat$algo <- gsub(" NA", "", x = dat$algo)
unique(dat$algo)

dat$algoF <- factor(dat$algo, levels = rev(c("rw", "stepping", "latent",
                                             "imh AUC_samples", "imh Laplace_analytic", "imh Laplace_analytic_wide",
                                             "gess AUC_samples", "gess Laplace_analytic", "gess Laplace_analytic_wide",
                                             "Qslice AUC_samples", "Qslice Laplace_analytic", "Qslice Laplace_analytic_wide")),
                    labels = rev(c("Random walk", "Steping out & shrinkage", "Latent slice",
                                   "Independence M-H: AUC", "Independence M-H: Laplace", "Independence M-H: Laplace (wide)",
                                   "Generalized elliptical slice: AUC", "Generalized elliptical slice: Laplace", "Generalized elliptical slice: Laplace (wide)",
                                   "Quantile slice: AUC", "Quantile slice: Laplace", "Quantile slice: Laplace (wide)"))
)

dat$n_extraF <- as.factor(dat$n_extra)

dat$typeF <- case_match(dat$type, c("gess", "latent", "stepping") ~ "Slice",
                        "imh" ~ "IMH", "Qslice" ~ "Quantile slice",
                        "rw" ~ "Rand walk") %>%
  factor(., levels = c("Rand walk", "Slice", "Quantile slice", "IMH"))

dat$targetlab <- case_match(dat$target, "hyper-g" ~ "h-g",
                            "hyper-g-log" ~ "h-g log")

dat$target_tx <- ifelse(grepl("log", dat$target), "Log transform", "Original") %>%
  as.factor()
dat$target_tx_alpha <- ifelse(dat$target_tx == "Log transform", 0.5, 1.0)


## evaluations per iteration
eval_tab <- dat %>% group_by(target, type, subtype, n_extra) %>% summarize(eval_per_iter_mean = mean(nEval / n_iter),
                                                      eval_per_iter_sd = sd(nEval / n_iter))

print(eval_tab, n = 30)

dat_plt <- dat %>% filter(n_extra %in% c(0))

## timing results
plt <- ggplot(dat_plt %>% filter(target %in% c("hyper-g")),
              aes(x = sampPsec / 1e3, y = algoF, fill = typeF), color = "gray") +
  geom_violin(draw_quantiles = 0.5 , scale = "width") +
  theme_bw() + theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  xlab("Effective samples (in thousands) per second") + ylab("") + labs(color = "", fill = "")

plt

ggsave(plot = plt, width = 6, height = 6,
       filename = paste0("plots/ESPS_primary_", dte, ".pdf"))


plt <- ggplot(dat_plt %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = sampPsec / 1e3, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  # geom_boxplot() +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  xlab("Effective samples (in thousands) per second") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 6, height = 6,
       filename = paste0("plots/ESPS_w_tnx_", dte, ".pdf"))



dat_plt_nxt <- dat %>% filter(n_extra %in% c(0, 100)) %>% 
  mutate(n_extra_lab = factor(n_extra, 
                              levels = c(0, 100), 
                              labels = c("Original Target", "Expensive Target")))

format_mean_sd <- function(mn, std) {
  # Round std for conditional logic
  std_fmt <- ifelse(std < 0.005 & std != 0, 0.01, std)
  
  # Choose format string
  ifelse(
    std == 0,
    sprintf("%.0f ± %.0f", mn, std),
    sprintf("%.1f ± %.2f", mn, std_fmt)
  )
}

stat_labels <- dat_plt_nxt %>%
  filter(target %in% c("hyper-g", "hyper-g-log")) %>%
  group_by(algoF, target_tx, n_extra_lab) %>%
  summarize(mEval = mean(nEval / n_iter), 
            sdEval = sd(nEval / n_iter),
            maxSampPsec = max(sampPsec),
            .groups = "drop") %>%
  mutate(algoF_num = as.numeric(factor(algoF, levels = sort(unique(dat_plt$algoF)))),
         dodge_offset = ifelse(target_tx == "Original", 0.2, -0.2),
         y = algoF_num + dodge_offset,
         x = Inf,
         label = format_mean_sd(mEval, sdEval)
         )

x_pad <- dat_plt_nxt %>% # hack to start the horizontal axis at 0
  group_by(n_extra_lab) %>%
  slice(1) %>%  # grab one row per facet
  mutate(sampPsec = 0)  # force x = 0

dat_plt_padded <- bind_rows(dat_plt_nxt, x_pad)

plt <- ggplot(dat_plt_nxt %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = sampPsec / 1.0, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + theme(legend.position = "none") +
  scale_alpha(guide = "none") +
  xlab("Effective samples per second") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black")) + 
  geom_text(data = stat_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, size = 3)

plt_nxt <- plt + facet_grid(cols = vars(n_extra_lab), scales = "free") + 
  coord_cartesian(clip = "off") + theme(panel.spacing = unit(3.5, "lines")) + 
  theme(plot.margin = unit(c(5.5, 50, 5.5, 0), "pt")) + 
  theme(panel.border = element_blank()) + geom_blank(data = x_pad, aes(x = sampPsec)) + 
  theme(strip.background = element_rect(fill = "gray90", color = NA))
plt_nxt

ggsave(plot = plt_nxt, width = 8, height = 4.5,
       filename = paste0("plots/ESPS_w_tnx_nextra", dte, ".pdf"))




## diagnostic

plt <- ggplot(dat_plt %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = g_mn, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Posterior mean: g") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/pm_g", dte, ".pdf"))

plt <- ggplot(dat_plt %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = g_sd, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Posterior SD: g") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

plt <- ggplot(dat_plt %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = g_SE, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") + xlim(c(0, 0.2)) +
  xlab("MCMC std. error: g") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/SE_g", dte, ".pdf"))

plt <- ggplot(dat_plt %>% filter(target %in% c("hyper-g", "hyper-g-log")),
              aes(x = nEval / n_iter, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Target evaluations per iteration") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

dat_plt %>% group_by(type, subtype) %>% summarize(mn_eval = mean(nEval/n_iter), sd_eval = sd(nEval/n_iter))

ggsave(plot = plt, width = 8, height = 8,
       filename = paste0("plots/nEval_", dte, ".pdf"))

plt <- ggplot(dat_plt %>% filter(target %in% c("hyper-g", "hyper-g-log"), type == "Qslice"),
              aes(x = auc, y = algoF, fill = typeF, color = target_tx, alpha = target_tx_alpha)) +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  theme_bw() + scale_alpha(guide = "none") +
  xlab("Transformed pseudo-target: AUC") + ylab("") + labs(color = "", fill = "") +
  scale_color_manual(values = c("gray", "black"))

plt

ggsave(plot = plt, width = 6, height = 5,
       filename = paste0("plots/AUC_", dte, ".pdf"))

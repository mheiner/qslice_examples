rm(list = ls())
dev.off()

set.seed(1)

library("qslice")
library("coda")
library("ggplot2")

source("0_data.R")
source("0_prior.R")
source("functions/full_conditionals.R")
source("functions/mcmc_gprior.R")
source("functions/tune.R")
source("functions/Laplace_g.R")
source("functions/invgamma.R")

state <- list(iter = 0)
state$psi <- prior$a_0 / prior$b_0
state$g <- prior$g_max / 2
state$beta <- prior$beta_0
state

n_iter <- 10e3

sampler <- list()

## plotting params
ylim_hist <- c(0.0, 2.1)
mar_use <- c(3.1, 3.1, 3.1, 0.4)
mgp_use <- c(1.5, 0.4, 0)
ylim_use <- c(0, 0.11)
height_use <- 2.75
width_use <- 4

### AUC Samples

## gather preliminary samples
sampler$g <- list(type = "stepping", subtype = NA, logG = F, w = 70.0)

mc_out <- mcmc_gprior(
  state = state,
  prior = prior,
  data = dat,
  sampler = sampler,
  n_iter = 2e3,
  save = TRUE,
  prog = 0
)

(state <- mc_out$state)

tmp_pseu <- pseudo_opt(
  samples = sapply(mc_out$sims, function(x) x$g),
  type = "samples",
  family = "t",
  df = c(1, 5),
  lb = 0.0,
  ub = prior$g_max,
  plot = TRUE,
  nbins = 20
)

sampler$g <- list(
  type = "Qslice",
  subtype = "AUC_samples",
  logG = sampler$g$logG,
  pseudo = tmp_pseu$pseudo,
  t = tmp_pseu$pseu$t
)

mc_out <- mcmc_gprior(
  state = state,
  prior = prior,
  data = dat,
  sampler = sampler,
  n_iter = n_iter,
  save = TRUE,
  prog = 0
)

(state <- mc_out$state)

draws_g <- sapply(mc_out$sims, function(x) x$g)
draws_u <- sapply(mc_out$extras, function(x) x$u)

gthm <- theme(
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 11)
)

ggplot(data.frame(u = draws_u), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  ylim(ylim_hist) +
  xlab(expression(psi)) +
  theme_classic() +
  gthm +
  labs(
    title = "AUC - samples",
    subtitle = bquote(
      "AUC: " ~ .(round(
        qslice::utility_shrinkslice(u = draws_u, utility_type = "AUC"),
        2
      )) ~
        "     " ~ IAT^{
        -1
      } ~ ": " ~ .(round(effectiveSize(draws_u) / n_iter, 3))
    )
  )
ggsave(file = "plots/histU_AUCsamples.pdf", height = 2.5, width = 3)

pdf(
  file = "plots/post_g_pseu_AUCsamples.pdf",
  height = height_use,
  width = width_use
)
par(mar = mar_use)
plot(
  density(draws_g),
  xlim = c(0, 70),
  xlab = expression(gamma),
  ylab = "density",
  main = "",
  bty = "n",
  mgp = mgp_use,
  ylim = ylim_use
)
curve(sampler$g$pseudo$d(x), from = 0, to = 75, add = TRUE, lty = 2, n = 1000)
legend(
  "topright",
  lty = c(1, 2),
  legend = c("Marginal posterior", "Pseudo-target"),
  bty = "n"
)
dev.off()

### Laplace Analytic

sampler$g <- list(
  type = "Qslice",
  subtype = "Laplace_analytic",
  sc_adj = 1.0,
  df = 1,
  logG = F
)

mc_out <- mcmc_gprior(
  state = state,
  prior = prior,
  data = dat,
  sampler = sampler,
  n_iter = n_iter,
  save = TRUE,
  prog = 0
)

(state <- mc_out$state)

draws_u <- sapply(mc_out$extras, function(x) x$u)

ggplot(data.frame(u = draws_u), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  ylim(ylim_hist) +
  xlab(expression(psi)) +
  theme_classic() +
  gthm +
  labs(
    title = "Laplace analytic",
    subtitle = bquote(
      "AUC: " ~ .(round(
        qslice::utility_shrinkslice(u = draws_u, utility_type = "AUC"),
        2
      )) ~
        "     " ~ IAT^{
        -1
      } ~ ": " ~ .(round(effectiveSize(draws_u) / n_iter, 3))
    )
  )
ggsave(file = "plots/histU_LaplaceAnalytic.pdf", height = 2.5, width = 3)

(beta_use <- sapply(mc_out$sims, function(x) x$beta) |> rowMeans())
(psi_use <- sapply(mc_out$sims, function(x) x$psi) |> mean())
qq0 <- dat$X %*% (beta_use - prior$beta_0) |> drop()
qq1 <- drop(crossprod(qq0))

tmp_pseu <- lapxt_g(
  p = dat$p,
  psi = psi_use,
  Q = qq1,
  a = prior$a_g,
  sc_adj = sampler$g$sc_adj,
  df = sampler$g$df,
  lb = 0.0,
  ub = prior$g_max,
  logG = sampler$g$logG
)

ltarget <- function(gg) {
  ## Likelihood for g is unnormalized Gaussian density for beta where
  ##  mean is beta_0 and
  ##  Cov is inv_XtX * g / psi
  ## Prior on g is hyper-g (Liang et al, 2008) on (0, g_max)

  if ((gg > 0.0) & (gg <= prior$g_max)) {
    qq <- psi_use * qq1 / gg
    logdet <- dat$p * log(gg) # only this part is a fn of g
    lpri <- prior$a_g * log(gg + 1.0)
    out <- -0.5 * (logdet + lpri + qq)
  } else {
    out <- -Inf
  }

  out
}

target <- function(gg) {
  exp(sapply(gg, FUN = ltarget))
}

int_tg <- integrate(target, lower = 0.0, upper = prior$g_max)
normc <- int_tg$value

gg_seq <- seq(0.0, 70.0, length = 1000)
fc_seq <- target(gg_seq) / normc

pdf(
  file = "plots/post_g_pseu_Laplace.pdf",
  height = height_use,
  width = width_use
)
par(mar = mar_use)
plot(
  gg_seq,
  fc_seq,
  type = "l",
  xlab = expression(gamma),
  ylab = "density",
  main = "",
  bty = "n",
  mgp = mgp_use,
  # ylim = c(0, max(fc_seq) * 1.1)
  ylim = ylim_use
)
curve(tmp_pseu$d(x), from = 0, to = 75, add = TRUE, lty = 2, n = 1000)
legend(
  "topright",
  lty = c(1, 2),
  legend = c("Full conditional", "Pseudo-target"),
  bty = "n"
)
dev.off()


### Laplace Analytic Wide

sampler$g <- list(
  type = "Qslice",
  subtype = "Laplace_analytic_wide",
  sc_adj = 1.5,
  df = 1,
  logG = F
)

mc_out <- mcmc_gprior(
  state = state,
  prior = prior,
  data = dat,
  sampler = sampler,
  n_iter = n_iter,
  save = TRUE,
  prog = 0
)

(state <- mc_out$state)

draws_u <- sapply(mc_out$extras, function(x) x$u)

ggplot(data.frame(u = draws_u), aes(x = u)) +
  geom_histogram(fill = "gray", aes(y = ..density..)) +
  ylim(ylim_hist) +
  xlab(expression(psi)) +
  theme_classic() +
  gthm +
  labs(
    title = "Laplace analytic - wide",
    subtitle = bquote(
      "AUC: " ~ .(round(
        qslice::utility_shrinkslice(u = draws_u, utility_type = "AUC"),
        2
      )) ~
        "     " ~ IAT^{
        -1
      } ~ ": " ~ .(round(effectiveSize(draws_u) / n_iter, 3))
    )
  )
ggsave(file = "plots/histU_LaplaceAnalyticWide.pdf", height = 2.5, width = 3)

## reuse beta and psi from previous
# (beta_use <- sapply(mc_out$sims, function(x) x$beta) |> rowMeans())
# (psi_use <- sapply(mc_out$sims, function(x) x$psi) |> mean())

qq0 <- dat$X %*% (beta_use - prior$beta_0) |> drop()
qq1 <- drop(crossprod(qq0))

tmp_pseu <- lapxt_g(
  p = dat$p,
  psi = psi_use,
  Q = qq1,
  a = prior$a_g,
  sc_adj = sampler$g$sc_adj,
  df = sampler$g$df,
  lb = 0.0,
  ub = prior$g_max,
  logG = sampler$g$logG
)

ltarget <- function(gg) {
  ## Likelihood for g is unnormalized Gaussian density for beta where
  ##  mean is beta_0 and
  ##  Cov is inv_XtX * g / psi
  ## Prior on g is hyper-g (Liang et al, 2008) on (0, g_max)

  if ((gg > 0.0) & (gg <= prior$g_max)) {
    qq <- psi_use * qq1 / gg
    logdet <- dat$p * log(gg) # only this part is a fn of g
    lpri <- prior$a_g * log(gg + 1.0)
    out <- -0.5 * (logdet + lpri + qq)
  } else {
    out <- -Inf
  }

  out
}

target <- function(gg) {
  exp(sapply(gg, FUN = ltarget))
}

(int_tg <- integrate(target, lower = 0.0, upper = prior$g_max))
normc <- int_tg$value

gg_seq <- seq(0.0, 70.0, length = 1000)
fc_seq <- target(gg_seq) / normc

pdf(
  file = "plots/post_g_pseu_LaplaceWide.pdf",
  height = height_use,
  width = width_use
)
par(mar = mar_use)
plot(
  gg_seq,
  fc_seq,
  type = "l",
  xlab = expression(gamma),
  ylab = "density",
  main = "",
  bty = "n",
  mgp = mgp_use,
  ylim = ylim_use
)
curve(tmp_pseu$d(x), from = 0, to = 75, add = TRUE, lty = 2, n = 1000)
legend(
  "topright",
  lty = c(1, 2),
  legend = c("Full conditional", "Pseudo-target"),
  bty = "n"
)
dev.off()

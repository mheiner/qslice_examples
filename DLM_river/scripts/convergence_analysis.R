
convergence_analysis <- function(fit, verbose = TRUE, multivariate = FALSE) {
  
  fit$mcmc_list <- lapply(fit$mcmc_list, function(x) as.mcmc(x)) |> mcmc.list()
  
  fit$mcmc_list$theta_dg <- lapply(fit$mcmc_list$theta, function(x) matrix(c(x), ncol = fit$parameters$n_iter) |> t() |> as.mcmc()) |> mcmc.list()

  if (exists("phi", fit$mcmc_list)) {
    fit$mcmc_list$phi_dg <- lapply(fit$mcmc_list$phi, function(x) matrix(c(x), ncol = fit$parameters$n_iter) |> t() |> as.mcmc()) |> mcmc.list()
  }
  
  fit$mcmc_list$sigma2_dg <- lapply(fit$mcmc_list$sigma2, function(x) matrix(x, ncol = 1) |> as.mcmc()) |> mcmc.list()
  fit$mcmc_list$phases_dg <- lapply(fit$mcmc_list$phases, function(x) t(x) |> as.mcmc()) |> mcmc.list()
  
  gelman_theta <- gelman.diag(fit$mcmc_list$theta_dg, autoburnin = FALSE, multivariate = multivariate)
  gelman_sig2 <- gelman.diag(fit$mcmc_list$sigma2_dg, autoburnin = FALSE, multivariate = multivariate)
  gelman_phase <- gelman.diag(fit$mcmc_list$phases_dg, autoburnin = FALSE, multivariate = multivariate)
  
  ess_theta <- effectiveSize(fit$mcmc_list$theta_dg)
  ess_sig2 <- effectiveSize(fit$mcmc_list$sigma2_dg)
  ess_phase <- effectiveSize(fit$mcmc_list$phases_dg)
  
  if (isTRUE(verbose)) {
    cat("Max Upper CI Gelman diagnostic, thetas", round(max(gelman_theta$psrf[,2]), 4), "\n")
    cat("Mean Upper CI Gelman diagnostic, thetas", round(mean(gelman_theta$psrf[,2]), 4), "\n")
    cat("Gelman diagnostic, phase\n")
    cat(round(gelman_phase$psrf[,1], 4), "\n")
    cat("Gelman diagnostic, phase upperCI\n")
    cat(round(gelman_phase$psrf[,2], 4), "\n")
    cat("Gelman diagnostic, err. variance\n")
    cat(round(gelman_sig2$psrf, 4), "\n")
    cat("Min effective sample size, theta:", min(ess_theta), "\n")
    cat("Avg effective sample size, theta:", mean(ess_theta), "\n")
  }
  
  list(gelman = list(theta = gelman_theta, sig2 = gelman_sig2, phase = gelman_phase),
       ess = list(theta = ess_theta, sig2 = ess_sig2, phase = ess_phase),
       mc_list = fit$mcmc_list)
}

QsliceFF_update <- function(ystar, f, m0, v0, tau2, sigma2, old_theta, time_diff, pseu_degf) {
  
  require("cucumber")
  
  # Put parameters into a list to be supplied to helper functions
  parameters <- list(
    ystar = ystar,
    f = f,
    m0 = m0,
    v0 = v0,
    tau2 = tau2,
    sigma2 = sigma2,
    old_theta = old_theta,
    time_diff = time_diff,
    pseu_degf = pseu_degf
  )
  
  n <- length(ystar)
  rejections <- 0
  
  kfilt <- ff(parameters = parameters) # this gives a sequence of moments based on old_theta
  stopifnot(length(kfilt$filtered_means) == n)

  ltarg <- function(x) {
    dlm_likelihood(theta = x, parameters = parameters)
  }
  
  ps_lpdf <- list()
  ps_qf <- list()
  ps_cdf <- list()
  
  for (i in 1:n) {
    tmp <- pseudo_t_list(loc = kfilt$filtered_means[i],
                         sc = sqrt(kfilt$filtered_variance[i]) * 1.0,
                         degf = parameters$pseu_degf,
                         lb = 0.0, ub = Inf)
    ps_lpdf[[i]] <- tmp$ld
    ps_qf[[i]] <- tmp$q
    ps_cdf[[i]] <- tmp$p
  }
  
  tmp <- slice_mv_transform(x = old_theta, 
                            target = ltarg, 
                            pseudo_log_pdf = ps_lpdf,
                            pseudo_inv_cdf = ps_qf,
                            pseudo_cdf = ps_cdf)

  list(theta = tmp$x, 
       rejections = tmp$nEvaluations - 2, 
       phi = tmp$u)
}

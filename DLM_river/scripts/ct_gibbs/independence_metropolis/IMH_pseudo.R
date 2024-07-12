imh_pseudo <- function(ystar, f, m0, v0, tau2, sigma2, old_theta, time_diff, pseu_family, pseu_params) {
  
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
    pseu_family = pseu_family, 
    pseu_params = pseu_params
  )
  
  # Extract length of the time series
  n <- length(ystar)
  
  # Calculate Kalman filter
  kfilt <- ff(parameters = parameters) # this gives a sequence of moments based on old_theta
  stopifnot(length(kfilt$filtered_means) == n)
  
  ltarg <- function(x) {
    dlm_likelihood(theta = x, parameters = parameters)
  }
  
  init_pseu_params <- pseu_params
  init_pseu_params$loc <- kfilt$filtered_means[n]
  init_pseu_params$sc <- sqrt(kfilt$filtered_variance[n]) * pseu_params$sc_adj
  
  pseu_cntl <- list( # this sequence will define the backward sampling parameters, so 
    
    pseu_init = pseudo_list(family = pseu_family,
                            params = init_pseu_params,
                            lb = 0.0, ub = Inf),
    
    loc_fn = function(z) { # function takes all z's up to (but not including) that point, i.e., z = (z[1], z[2]) = (theta[n], theta[n-1])
      k <- length(z) + 1
      tt <- n - k + 1
      theta_next <- z[k - 1] # most recent in the vector passed in
      out <- kfilt$filtered_means[tt] + (kfilt$filtered_variance[tt] * (theta_next - kfilt$prior_mean[tt+1]) / kfilt$prior_variance[tt+1])
      out
    },
    
    sc_fn = function(z) {
      k <- length(z) + 1
      tt <- n - k + 1
      var_out <- kfilt$filtered_variance[tt] * (1.0 - kfilt$filtered_variance[tt] / kfilt$prior_variance[tt+1])
      sqrt(var_out) * pseu_params$sc_adj
    },
    
    lb = rep(0.0, n),
    ub = rep(Inf, n)
  )
  
  tmp_seq <- pseudo_condseq(x = rev(old_theta), 
                            pseudo_init = pseu_cntl$pseu_init,
                            loc_fn = pseu_cntl$loc_fn,
                            sc_fn = pseu_cntl$sc_fn,
                            lb = pseu_cntl$lb,
                            ub = pseu_cntl$ub)
  
  lfx0 <- ltarg(old_theta) - sum(sapply(1:n, function(i) tmp_seq[[i]]$ld(old_theta[n-i+1])))
  
  u1 <- runif(n, min = 0.0, max = 1.0)
  tmp <- pseudo_condseq_XfromU(u = u1, 
                               pseudo_init = pseu_cntl$pseu_init,
                               loc_fn = pseu_cntl$loc_fn,
                               sc_fn = pseu_cntl$sc_fn,
                               lb = pseu_cntl$lb,
                               ub = pseu_cntl$ub)
  
  theta_proposal <- tmp$x |> rev()
  tmp_seq1 <- tmp$pseudo_seq
  
  lfx1 <- ltarg(theta_proposal) - sum(sapply(1:n, function(i) tmp_seq1[[i]]$ld(tmp$x[i])))
  
  # Calculate log acceptance rate
  log_acceptance <- lfx1 - lfx0
  
  # Log cutoff draw
  z <- log(runif(1))
  
  # Accept or reject
  accept <- log_acceptance > z
  
  # Return the accepted value or the old value
  if(accept) {
    return(theta_proposal)
  } else {
    return(old_theta)
  }
}

### Doesn't really work--increases rejections and is slower

Qslice_targetShrink_update <- function(ystar, f, m0, v0, tau2, sigma2, old_theta, time_diff, pseu_degf) {
  
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

  ltarg <- function(x) { # slice sampler will consider theta vector in reverse order; this function receives it in natural order
    llik_compon <- dlm_likelihood_components(theta = rev(x), parameters = parameters)
    out <- list()
    out$components <- llik_compon$llik_proc + llik_compon$llik_obs
    out
  }
  
  pseu_cntl <- list( # this sequence will define the backward sampling parameters, so 
    
    pseu_init = pseudo_t_list(loc = kfilt$filtered_means[n], 
                              sc = sqrt(kfilt$filtered_variance[n]), 
                              degf = pseu_degf,
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
      sqrt(var_out)
    },
    
    degf = pseu_degf,
    lb = rep(0.0, n),
    ub = rep(Inf, n)
  )
  
  tmp <- slice_mv_transform_seq_targetShrink(x = rev(old_theta), 
                                target = ltarg, 
                                pseudo_control = pseu_cntl)

  list(theta = rev(tmp$x), 
       rejections = tmp$nEvaluations - 2, 
       phi = rev(tmp$u))
}
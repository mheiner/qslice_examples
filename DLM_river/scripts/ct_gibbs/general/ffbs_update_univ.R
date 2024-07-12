ffbs_update <- function(ystar, f, m0, v0, tau2, sigma2, time_diff) {
  
  kfilt <- ff(ystar = ystar, f = f, m0 = m0, v0 = v0, 
              tau2 = tau2, sigma2 = sigma2,
              time_diff = time_diff)
  
  backward_sample(kfilt)
}

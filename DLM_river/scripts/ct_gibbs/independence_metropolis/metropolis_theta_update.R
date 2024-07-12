update_theta <- function() {
  
  for(i in 1:p) {
    
    # Use the current thetas and F to create a ystar time series
    ystar <- y - colSums(t(F_mat[,-i]) * theta_now[-i,, drop = FALSE])
    
    # Update theta_i with the appropriate method
    if(i != 1){
      # theta_now[i,] <- imh_pseudo(ystar, F_mat[, i], m0[i], C0[i, i], W[i, i],
      #                             sigma2_now, theta_now[i,], time_diff,
      #                             pseu_family, pseu_params)
      theta_now[i,] <- independence_metropolis_update(ystar, F_mat[, i], m0[i], C0[i, i], W[i, i],
                                  sigma2_now, theta_now[i,], time_diff,
                                  pseu_family, pseu_params)
    } else {
      theta_now[i,] <- ffbs_update(ystar, F_mat[, i], m0[i], C0[i, i], W[i, i], sigma2_now, time_diff)
    }
  }
  
  list(theta = theta_now, rejections = 0) # these refer to rejections that require another sample
}

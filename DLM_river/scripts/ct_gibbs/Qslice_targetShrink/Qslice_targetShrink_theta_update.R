update_theta <- function() {
  
  rejections <- 0
  phi_mat <- matrix(NA, nrow = n_freq, ncol = n)
  
  for (i in 1:p) {
    
    # Use the current thetas and F to create a ystar time series
    ystar <- y - colSums(t(F_mat[,-i]) * theta_now[-i,, drop = FALSE])
    
    # Update theta_i with the appropriate method
    if (i != 1) {
      tmp <- Qslice_targetShrink_update(ystar, F_mat[, i], m0[i], C0[i, i], W[i, i], sigma2_now, theta_now[i,], time_diff, pseu_degf)
      theta_now[i,] <- tmp$theta
      rejections <- rejections + tmp$rejections
      phi_mat[i-1,] <- tmp$phi
    } else {
      theta_now[i,] <- ffbs_update(ystar, F_mat[, i], m0[i], C0[i, i], W[i, i], sigma2_now, time_diff)
    }
  }
  
  list(theta = theta_now, rejections = rejections, phi_mat = phi_mat)
}

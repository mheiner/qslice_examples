update_theta <- function() {
  
  rejections <- 0

  for(i in 1:p) {
    
    # Use the current thetas and F to create a ystar time series
    ystar <- y - colSums(t(F_mat[,-i]) * theta_now[-i,, drop = FALSE])
    
    # Update theta_i with the appropriate method
    if(i != 1) {
      latent_slice <- latent_slice_update(ystar, F_mat[, i], m0[i], C0[i, i], W[i, i], 
                                          sigma2_now, s_mat[i-1,], latent_scale, theta_now[i,], time_diff)
      rejections <- rejections + latent_slice$rejections
      theta_now[i,] <- latent_slice$theta
      s_mat[i-1,] <- latent_slice$s
    } else {
      theta_now[i,] <- ffbs_update(ystar, F_mat[, i], m0[i], C0[i, i], W[i, i], sigma2_now, time_diff)
    }
    
  }
  
  list(theta = theta_now, rejections = rejections, s_mat = s_mat)
}

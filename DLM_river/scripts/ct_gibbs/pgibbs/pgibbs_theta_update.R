update_theta <- function() {

  for (i in 1:p) {

    # Use the current thetas and F to create a ystar time series
    ystar <- y - colSums(t(F_mat[,-i]) * theta_now[-i,, drop = FALSE])

    # Update theta_i with the appropriate method
    if(i != 1){
      theta_now[i,] <- particle_filter_t1_update(n_particles = n_particles, 
                                       ystar = ystar, 
                                       f = F_mat[, i], 
                                       m0 = m0[i], 
                                       v0 = C0[i, i], 
                                       tau2 = W[i, i], 
                                       sigma2 = sigma2_now, 
                                       old_theta = theta_now[i,], 
                                       time_diff = time_diff)
      
    } else {
      theta_now[i,] <- ffbs_update(ystar, F_mat[, i], m0[i], C0[i, i], W[i, i], sigma2_now, time_diff)
    }
  }
  
  list(theta = theta_now, rejections = 0)
}

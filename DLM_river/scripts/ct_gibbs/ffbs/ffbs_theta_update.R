update_theta <- function() {
  theta_new <- ffbs_rw(y = y, F_mat = F_mat, m0 = m0, C0 = C0, W = W, 
                       sigma2 = sigma2_now, time_diff = time_diff)
  list(theta = theta_new, rejections = 0)
}

### Create function to update variance

sigma2_update <- function(y, a, b, F_mat, theta_matrix) {
  n <- length(y)
  yhat <- colSums(t(F_mat) * theta_matrix)
  
  astar <- a + n/2
  bstar <- b + 0.5 * crossprod(y - yhat)
  
  full_conditional_draw <- 1.0 / rgamma(1, shape = astar, rate = bstar)
  
  return(full_conditional_draw)
}

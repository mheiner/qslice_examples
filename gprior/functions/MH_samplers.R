random_walk_sampler <- function(lf, support, x_0, c) {
  x_out <- x_0
  accept <- 0
  n_eval <- 0

  ## proposed draw
  x.dot <- rnorm(1, x_0, c)

  if (x.dot >= support[1] && x.dot <= support[2]) {
    logr <- lf(x.dot) - lf(x_0)
    u <- runif(1, 0, 1)
    n_eval <- n_eval + 2

    if (log(u) < logr) {
      x_out <- x.dot
      accept <- 1
    }
  }

  list(x = x_out, accept = accept, nEvaluations = n_eval)
}

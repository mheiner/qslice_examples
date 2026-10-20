doextra <- function(n) {
  Z <- matrix(rnorm(2 * n^2), ncol = n)
  ZtZ <- crossprod(Z)
  CZtZ <- chol(ZtZ)
  chol2inv(CZtZ)
}

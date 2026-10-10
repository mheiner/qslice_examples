lCC <- matrix(rnorm(n * K, mean = C_mu, sd = C_sig), nrow = n)
CC <- exp(lCC)
lCCsum <- apply(lCC, 1L, logsumexp)
CCsum <- exp(lCCsum)

theta <- rgamma(K, shape = a0)
(theta <- theta / sum(theta))

if (type == "decreasing") {
  theta <- sort(theta, decreasing = TRUE)
}

thetaC <- rep(theta, each = n) * CC

x <- numeric(n)
for (i in 1:n) {
  x[i] <- sample.int(K, size = 1, prob = thetaC[i, ])
}
x

(xtab <- tabulate(x, nbins = K))
xtab / n

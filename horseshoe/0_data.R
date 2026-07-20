dat <- list()

library("lars")
data(diabetes)

if (data_use == "db40") {
  set.seed(1)
  indx_obs <- sort(sample.int(length(diabetes$y), replace = FALSE, size = 40))
} else {
  indx_obs <- 1:length(diabetes$y)
}

dat$y <- scale(diabetes$y[indx_obs])
dat$X <- scale(diabetes$x2[indx_obs, ])

rm(diabetes)

dat$XtX <- crossprod(dat$X)

library("glmnet")
lso <- cv.glmnet(dat$X, dat$y, alpha = 1, family = "gaussian")
dat$beta_hat <- as.numeric(coef(lso)[-1])
str(dat$beta_hat)

(dat$n <- length(dat$y))
(dat$p <- ncol(dat$X))

dat

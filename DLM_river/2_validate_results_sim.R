
## To be run on MCMC outputs: fit and mc_diag

plot(mc_diag$mc_list$phases_dg)

plot(mc_diag$mc_list$sigma2_dg)

tt <- 10
k <- 2
thin_seq <- seq(1, n_iter, length = 5000) |> floor()
plot((fit$mcmc_list$theta[[1]][k,tt,thin_seq]), type = "l", lty = 2)
lines((fit$mcmc_list$theta[[2]][k,tt,thin_seq]), col = "red", lty = 2)


source("scripts/dhr_plotting.R")
fit <- combine_chains(fit)
str(fit)

dg <- diagnostics_dhr(y = fit$parameters$y, dhr_object = fit)


plot_fit <- plot_dhr(y = fit$parameters$y, dhr_object = fit)
plot_fit <- plot_fit  + 
  scale_y_continuous(sec.axis = sec_axis(~ exp(.), 
                                         name = TeX("$\\NO_3^-$  (mg / L)"), 
                                         breaks = c(1, 2, 5, 10, 20, 50)),
                     ) +
  geom_hline(yintercept = log(19), color = "gray50", lty = 2)  


plot(fit$parameters$times, fit$parameters$y, type = "l")
points(fit$parameters$times, dg$fitted_values, col = "red")
lines(fit$parameters$times, dg$fitted_values, col = "red")

acf(dg$residuals)
pacf(dg$residuals)


source("scripts/dlm_plotting.R")
plt_theta <- plot_theta(fit$mcmc_sample$theta_sample, times = fit$parameters$times,
                        theta_labels = c("level", "amplitude: annual", "amplitude: semiannual", "amplitude: triannual"),
                        theta_symbols = c("$\\beta_0$", "$\\alpha_1$", "$\\alpha_2$", "$\\alpha_3$"))
print(plt_theta)

library("gridExtra")

pdf("plotfit_river_long_pits_normal.pdf", width = 7.5, height = 4.75, onefile = TRUE)
grid.arrange(plot_fit, plt_theta[[1]], plt_theta[[2]], plt_theta[[3]],
             layout_matrix = rbind(c(1,1,1), c(2,3,4)))
dev.off()

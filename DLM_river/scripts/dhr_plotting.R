####################
### DHR Plotting ###
####################

library(latex2exp)
library(magrittr)
library(ggplot2)


combine_chains <- function(dhr_fit) {
  
  ### Shapes we need from original code
  # mcmc_sample <- list(
  #   theta_sample = array(dim = c(p, n, n_iter * n_chains)),
  #   sigma2_sample = numeric(n_iter * n_chains),
  #   phases_sample = matrix(nrow = n_freq, ncol = n_iter * n_chains)
  # )
  
  
  dhr_fit$mcmc_sample <- list()

  ## matrices with rows as iterations  
  # tmp <- lapply(dhr_fit$mcmc_list$theta, function(x) matrix(c(x), ncol = fit$parameters$n_iter) |> t())
  # dhr_fit$mcmc_sample$theta_sample <- do.call(rbind, tmp)
  
  dhr_fit$mcmc_sample$theta_sample <- array(do.call(c, dhr_fit$mcmc_list$theta), 
                                            c(length(fit$parameters$m0), 
                                              length(fit$parameters$times), 
                                              fit$parameters$n_chains*fit$parameters$n_iter))
  
  dhr_fit$mcmc_sample$phases_sample <- do.call(cbind, dhr_fit$mcmc_list$phases)

  # tmp <- lapply(dhr_fit$mcmc_list$sigma2, function(x) matrix(x, ncol = 1)) 
  dhr_fit$mcmc_sample$sigma2_sample <- do.call(c, dhr_fit$mcmc_list$sigma2)
  
  dhr_fit
}


diagnostics_dhr <- function(y, dhr_object, state = NULL, conf_level = 0.9) {
  
  # Pull the period and frequencies from the dhr_object
  period <- dhr_object$parameters$period
  frequencies <- dhr_object$parameters$frequencies
  
  # Pull the thetas and phases from the dhr_object
  theta_sample <- dhr_object$mcmc_sample$theta_sample
  phases_sample <- dhr_object$mcmc_sample$phases_sample
  
  # Check if the dhr_object is continuous time
  continuous_time <- "times" %in% names(dhr_object$parameters)
  if(continuous_time) {
    times <- dhr_object$parameters$times
  }
  
  ### Calculate mean and bounds based on draws ###
  
  # Extract iterations
  J <- dim(theta_sample)[3]
  # Calculate average theta from the drawn theta vectors
  posterior_theta_mean <- t(apply(theta_sample, MARGIN = c(1,2), FUN = mean))
  # Calculate number of observations
  n <- nrow(posterior_theta_mean)
  p <- ncol(posterior_theta_mean)
  # Find individual fitted observation values by converting
  # each sampled vector and phases to the observation scale
  posterior_fitted_values <- matrix(nrow = J, ncol = n)
  for(j in 1:J) {
    F_mat <- matrix(nrow = n, ncol = p)
    F_mat[,1] <- 1
    if(continuous_time){
      for(i in 2:p) {
        F_mat[,i] <- cos(2 * pi * (times * frequencies[i-1] + phases_sample[i-1, j]) / period)
      }
    } else {
      for(i in 2:p) {
        F_mat[,i] <- cos(2 * pi * ((1:n) * frequencies[i-1] + phases_sample[i-1, j]) / period)
      }
    }
    for(i in 1:n) {
      posterior_fitted_values[j,i] <- crossprod(F_mat[i,], theta_sample[,i,j])
    }
  }
  
  posterior_fitted_mean <- colMeans(posterior_fitted_values)
  
  # Calculate quantiles from desired confidence level
  u_quant <- (1 + conf_level)/2
  l_quant <- (1 - conf_level)/2
  
  # Bounds based on sampling results using quantiles
  u_posterior_fitted_values <- apply(posterior_fitted_values, MARGIN = 2, FUN = quantile, u_quant)
  l_posterior_fitted_values <- apply(posterior_fitted_values, MARGIN = 2, FUN = quantile, l_quant)
  
  # Calculate values of interest
  if(!is.null(state)){
    coverage <- mean(state < u_posterior_fitted_values & state > l_posterior_fitted_values)
  }
  
  width <- mean(u_posterior_fitted_values - l_posterior_fitted_values)
  rmse <- (y - posterior_fitted_mean)^2 %>% mean() %>% sqrt()
  bias <- mean(y - posterior_fitted_mean)
  
  if(is.null(state)){
    diagnostics <- c(width = width,
                     rmse = rmse,
                     bias = bias)
  } else {
    diagnostics <- c(coverage = coverage,
                     width = width,
                     rmse = rmse,
                     bias = bias)
  }
  output <- list(residuals = y - posterior_fitted_mean,
                 diagnostics = diagnostics,
                 fitted_values = posterior_fitted_mean)
  
  print(diagnostics)
  
  return(output)
}



plot_dhr <- function(y = NULL, dhr_object, state = NULL, conf_level = 0.9) {
  
  # Pull the period and frequencies from the dhr_object
  period <- dhr_object$parameters$period
  frequencies <- dhr_object$parameters$frequencies
  
  # Pull the thetas and phases from the dhr_object
  theta_sample <- dhr_object$mcmc_sample$theta_sample
  phases_sample <- dhr_object$mcmc_sample$phases_sample
  
  # Check if the dhr_object is continuous time
  continuous_time <- "times" %in% names(dhr_object$parameters)
  if(continuous_time) {
    times <- dhr_object$parameters$times
  }
  
  ### Calculate mean and bounds based on draws ###
  
  # Extract iterations
  J <- dim(theta_sample)[3]
  # Calculate average theta from the drawn theta vectors
  posterior_theta_mean <- t(apply(theta_sample, MARGIN = c(1,2), FUN = mean))
  # Calculate number of observations
  n <- nrow(posterior_theta_mean)
  p <- ncol(posterior_theta_mean)
  # Find individual fitted observation values by converting
  # each sampled vector and phases to the observation scale
  posterior_fitted_values <- matrix(nrow = J, ncol = n)
  for(j in 1:J) {
    F_mat <- matrix(nrow = n, ncol = p)
    F_mat[,1] <- 1
    if(continuous_time){
      for(i in 2:p) {
        F_mat[,i] <- cos(2 * pi * (times * frequencies[i-1] + phases_sample[i-1, j]) / period)
      }
    } else {
      for(i in 2:p) {
        F_mat[,i] <- cos(2 * pi * ((1:n) * frequencies[i-1] + phases_sample[i-1, j]) / period)
      }
    }
    for(i in 1:n) {
      posterior_fitted_values[j,i] <- crossprod(F_mat[i,], theta_sample[,i,j])
    }
  }
  
  posterior_fitted_mean <- colMeans(posterior_fitted_values)
  
  # Calculate quantiles from desired confidence level
  u_quant <- (1 + conf_level)/2
  l_quant <- (1 - conf_level)/2
  
  # Bounds based on sampling results using quantiles
  u_posterior_fitted_values <- apply(posterior_fitted_values, MARGIN = 2, FUN = quantile, u_quant)
  l_posterior_fitted_values <- apply(posterior_fitted_values, MARGIN = 2, FUN = quantile, l_quant)
  
  ### Now plot the smoothed posterior with error bars ###
  if(is.null(state)) {
    # Basic frame to plot time against responses
    plot_frame1 <- data.frame(tt = rep(1:n,2),
                              response = c(y, posterior_fitted_mean),
                              type = factor(c(rep("Data", n), rep("Fit",n)), levels = c("Fit", "Data")) )   
  } else {
    # Basic frame to plot time against responses
    plot_frame1 <- data.frame(tt = rep(1:n,2),
                              response = c(state, posterior_fitted_mean),
                              type = c(rep("theta",n),rep("posterior",n)))
  }
  
  # Frame for posterior error bands
  plot_frame2 <- data.frame(tt = c(1:n,rev(1:n)),
                            bounds = c(l_posterior_fitted_values, rev(u_posterior_fitted_values)))
  
    # Overwrite discrete time if we have continuous time data
  if(continuous_time) {
    plot_frame1$tt <- rep(times, 2)
    plot_frame2$tt <- c(times, rev(times))
  }
  
  # Initialize plot
  p <- ggplot() + theme_minimal()
  # Plot connecting lines color-coded by response type
  p <- p + geom_line(data = plot_frame1, mapping = aes(x = tt, y = response, color = type, linetype = type), linewidth = 1.0) +
    theme(legend.title = element_blank()) + 
    geom_point(plot_frame1 %>% filter(type == "Data"), mapping = aes(x = tt, y = response, color = type), size = 1.9)
  # Plot error bands
  p <- p + geom_polygon(data = plot_frame2, aes(x = tt, y = bounds), fill = "#F8766D", alpha = 0.3)
  # Clarify legend labels
  if(is.null(state)){
    # p <- p + scale_color_discrete(labels = c("Fit", "Data"))
    p <- p + theme(legend.position="none")
  } else {
    p <- p + scale_color_discrete(labels = c(TeX("$E(f\\theta|y)$"), TeX("True $f\\theta$")))
  }
  # Clarify plot axis and legend titles
  p <- p + labs(x = "Year", y = TeX("$\\log\\ NO_3^-$"), color = NULL) + 
    scale_y_continuous(sec.axis = sec_axis(~ exp(.), name = TeX("$\\NO_3^-$ in mg/L")))
  print(p)
}

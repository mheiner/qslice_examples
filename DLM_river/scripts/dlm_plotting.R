####################
### DLM Plotting ###
####################

library(latex2exp)
library(magrittr)
library(ggplot2)



plot_theta <- function(theta_sample, theta_index = NULL, theta = NULL, times = NULL, 
                       conf_level = .9, theta_labels = NULL, theta_symbols = NULL) {
  
  # Show all theta if no index is specified
  if(is.null(theta_index)){
    theta_index <- 1:nrow(theta_sample)
  }
  
  # Calculate useful quantities
  # Calculate quantiles from desired confidence level
  u_quant <- (1 + conf_level)/2
  l_quant <- (1 - conf_level)/2
  
  plots <- list()
  
  createplot <- function(theta_label = NULL, theta_symbol = NULL) {
    # Initialize plot
    p <- ggplot() + theme_minimal()
    # Plot connecting lines color-coded by response type
    p <- p + geom_line(data = plot_frame1, 
                       mapping = aes(x = tt, y = response, color = type), 
                       linewidth = 1.0)
    # Plot error bands
    p <- p + geom_polygon(data = plot_frame2, aes(x = tt, y = bounds), fill = "#F8766D", alpha = 0.3)
    # Clarify legend labels
    
    if (is.null(theta)) {
      p <- p + scale_color_discrete(labels = TeX("$E(\\theta|y)$"))
    } else {
      p <- p + scale_color_discrete(labels = c(TeX(paste("$E(\\theta_",index,"|y)$")), TeX(paste("True $\\theta_",index,"$"))))
    }
    
    # Clarify plot axis and legend titles
    p <- p + 
      labs(x = "", y = TeX(theta_symbol), color = "Type") + 
      theme(legend.position="none") + ggtitle(label = theta_label) + 
      theme(plot.title = element_text(size = 11))
    p
  }

  # Determine class of theta sample to see if we have a univariate or multivariate DLM
  theta_sample_class <- class(theta_sample)
  
  if (all(theta_sample_class == "array")){
    
    n <- ncol(theta_sample)
    
    for (index in theta_index) {
      # Subset to the current theta of interest
      current_theta <- theta_sample[index,,]
      
      # Calculate the mean
      current_theta_mean <- apply(current_theta, 1, FUN = mean)
      
      # Bounds based on sampling results using quantiles
      u_current_theta <- apply(current_theta, MARGIN = 1, FUN = quantile, u_quant)
      l_current_theta <- apply(current_theta, MARGIN = 1, FUN = quantile, l_quant)
      
      ### Now plot the smoothed posterior with error bars ###
      if (is.null(theta)) {
        # Basic frame to plot time against responses
        plot_frame1 <- data.frame(tt = 1:n,
                                  response = current_theta_mean,
                                  type = rep("posterior",n))
        # Frame for posterior credible bands
        plot_frame2 <- data.frame(tt = c(1:n,rev(1:n)),
                                  bounds = c(l_current_theta,
                                             rev(u_current_theta)))
        
        if(!is.null(times)) {
          plot_frame1$tt <- times
          plot_frame2$tt <- c(times, rev(times))
        }
        
      } else {
        # Extract actual values of theta
        current_true_theta <- theta[,index]
        # Basic frame to plot time against responses
        plot_frame1 <- data.frame(tt = rep(1:n,2),
                                  response = c(current_true_theta, current_theta_mean),
                                  type = c(rep("true",n), rep("posterior",n)))
        # Frame for posterior credible bands
        plot_frame2 <- data.frame(tt = c(1:n,rev(1:n)),
                                  bounds = c(l_current_theta,
                                             rev(u_current_theta)))
        
      }
      
      plots[[index]] <- createplot(theta_label = theta_labels[index], theta_symbol = theta_symbols[index])
    }
    
  } else if (any(theta_sample_class == "matrix")) {
    n <- nrow(theta_sample)
    
    # Calculate the mean
    theta_mean <- apply(theta_sample, 1, FUN = mean)
    
    # Bounds based on sampling results using quantiles
    u_theta <- apply(theta_sample, MARGIN = 1, FUN = quantile, u_quant)
    l_theta <- apply(theta_sample, MARGIN = 1, FUN = quantile, l_quant)
    
    ### Now plot the smoothed posterior with error bars ###
    if (is.null(theta)) {
      # Basic frame to plot time against responses
      plot_frame1 <- data.frame(tt = 1:n,
                                response = theta_mean,
                                type = rep("posterior",n))
      # Frame for posterior credible bands
      plot_frame2 <- data.frame(tt = c(1:n,rev(1:n)),
                                bounds = c(l_theta,
                                           rev(u_theta)))
      

    } else {
      # Basic frame to plot time against responses
      plot_frame1 <- data.frame(tt = rep(1:n,2),
                                response = c(theta, theta_mean),
                                type = c(rep("true",n), rep("posterior",n)))
      # Frame for posterior credible bands
      plot_frame2 <- data.frame(tt = c(1:n,rev(1:n)),
                                bounds = c(l_theta,
                                           rev(u_theta)))
      
    }
    
    if(!is.null(times)) {
      plot_frame1$tt <- times
      plot_frame2$tt <- c(times, rev(times))
    }
  
    plots <- createplot()
  }
  
  plots
}








diagnostics_DLM <- function(y, theta_sample, F_mat, state = NULL, conf_level = 0.9) {
  ### Calculate mean and bounds based on draws ###
  if(!is.matrix(theta_sample)) {
    # Extract iterations
    J <- dim(theta_sample)[3]
    # Calculate average theta from the drawn theta vectors
    posterior_theta_mean <- t(apply(theta_sample, MARGIN = c(1,2), FUN = mean))
    # Calculate number of observations
    n <- nrow(posterior_theta_mean)
    # Convert posterior mean vectors to posterior mean observations
    posterior_mean_obs <- numeric(n)
    for(i in 1:n){
      posterior_mean_obs[i] <- crossprod(F_mat[i,], posterior_theta_mean[i,])
    }
    # For bounds also find individual sampled observation values by converting
    # each sampled vector, not just the means
    posterior_obs <- matrix(nrow = J, ncol = n)
    for(i in 1:n) {
      for(j in 1:J) {
        posterior_obs[j,i] <- crossprod(F_mat[i,], theta_sample[,i,j])
      }
    }
  } else {
    # Extract iterations
    J <- ncol(theta_sample)
    # Calculate average theta from the drawn thetas
    posterior_theta_mean <- rowMeans(theta_sample)
    # Calculate number of observations
    n <- length(posterior_theta_mean)
    # Convert posterior mean vectors to posterior mean observations
    posterior_mean_obs <- numeric(n)
    for(i in 1:n){
      posterior_mean_obs[i] <- crossprod(F_mat[i,], posterior_theta_mean[i])
    }
    # For bounds also find individual sampled observation values by converting
    # each sampled vector, not just the means
    posterior_obs <- matrix(nrow = J, ncol = n)
    for(i in 1:n) {
      for(j in 1:J) {
        posterior_obs[j,i] <- crossprod(F_mat[i,], theta_sample[i,j])
      }
    }
  }
  
  # Calculate quantiles from desired confidence level
  u_quant <- (1 + conf_level)/2
  l_quant <- (1 - conf_level)/2
  
  # Bounds based on sampling results using quantiles
  u_posterior_obs <- apply(posterior_obs, MARGIN = 2, FUN = quantile, u_quant)
  l_posterior_obs <- apply(posterior_obs, MARGIN = 2, FUN = quantile, l_quant)
  
  # Calculate values of interest
  if(!is.null(state)){
    coverage <- mean(state < u_posterior_obs & state > l_posterior_obs)
  }
  width <- mean(u_posterior_obs - l_posterior_obs)
  rmse <- (y - posterior_mean_obs)^2 %>% mean() %>% sqrt()
  bias <- mean(y - posterior_mean_obs)
  
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
  output <- list(residuals = y - posterior_mean_obs,
                 diagnostics = diagnostics,
                 fitted_values = posterior_mean_obs)
  return(output)
}

plot_DLM <- function(y = NULL, theta_sample, F_mat, state = NULL, conf_level = 0.9) {
  ### Calculate mean and bounds based on draws ###
  
  if(!is.matrix(theta_sample)) {
    # Extract iterations
    J <- dim(theta_sample)[3]
    # Calculate average theta from the drawn theta vectors
    posterior_theta_mean <- t(apply(theta_sample, MARGIN = c(1,2), FUN = mean))
    # Calculate number of observations
    n <- nrow(posterior_theta_mean)
    # Convert posterior mean vectors to posterior mean observations
    posterior_mean_obs <- numeric(n)
    for(i in 1:n){
      posterior_mean_obs[i] <- crossprod(F_mat[i,], posterior_theta_mean[i,])
    }
    # For bounds also find individual sampled observation values by converting
    # each sampled vector, not just the means
    posterior_obs <- matrix(nrow = J, ncol = n)
    for(i in 1:n) {
      for(j in 1:J) {
        posterior_obs[j,i] <- crossprod(F_mat[i,], theta_sample[,i,j])
      }
    }
  } else {
    # Extract iterations
    J <- ncol(theta_sample)
    # Calculate average theta from the drawn thetas
    posterior_theta_mean <- rowMeans(theta_sample)
    # Calculate number of observations
    n <- length(posterior_theta_mean)
    # Convert posterior mean vectors to posterior mean observations
    posterior_mean_obs <- numeric(n)
    for(i in 1:n){
      posterior_mean_obs[i] <- crossprod(F_mat[i,], posterior_theta_mean[i])
    }
    # For bounds also find individual sampled observation values by converting
    # each sampled vector, not just the means
    posterior_obs <- matrix(nrow = J, ncol = n)
    for(i in 1:n) {
      for(j in 1:J) {
        posterior_obs[j,i] <- crossprod(F_mat[i,], theta_sample[i,j])
      }
    }
  }
  
  # Calculate quantiles from desired confidence level
  u_quant <- (1 + conf_level)/2
  l_quant <- (1 - conf_level)/2
  
  # Bounds based on sampling results using quantiles
  u_posterior_obs <- apply(posterior_obs, MARGIN = 2, FUN = quantile, u_quant)
  l_posterior_obs <- apply(posterior_obs, MARGIN = 2, FUN = quantile, l_quant)
  
  ### Now plot the smoothed posterior with error bars ###
  if(is.null(state)) {
    # Basic frame to plot time against responses
    plot_frame1 <- data.frame(tt = rep(1:n,2),
                              response = c(y,posterior_mean_obs),
                              type = c(rep("y",n),rep("posterior",n)))    
  } else {
    # Basic frame to plot time against responses
    plot_frame1 <- data.frame(tt = rep(1:n,2),
                              response = c(state,posterior_mean_obs),
                              type = c(rep("theta",n),rep("posterior",n)))
  }
  # Frame for posterior error bands
  plot_frame2 <- data.frame(tt = c(1:n,rev(1:n)),
                            bounds = c(l_posterior_obs,rev(u_posterior_obs)))
  
  # Initialize plot
  p <- ggplot() + theme_minimal()
  # Plot connecting lines color-coded by response type
  p <- p + geom_line(data = plot_frame1, mapping = aes(x = tt, y = response, color = type))
  # Plot error bands
  p <- p + geom_polygon(data = plot_frame2, aes(x = tt, y = bounds), fill = "#F8766D", alpha = 0.3)
  # Clarify legend labels
  if(is.null(state)){
    p <- p + scale_color_discrete(labels = c(TeX("$E(f\\theta|y)$"),"Observed"))
  } else {
    p <- p + scale_color_discrete(labels = c(TeX("$E(f\\theta|y)$"),TeX("True $f\\theta$")))
  }
  # Clarify plot axis and legend titles
  p <- p + labs(x = "Time", y = TeX("$\\log(NO_3^-)$"), color = "Type")
  print(p)
}

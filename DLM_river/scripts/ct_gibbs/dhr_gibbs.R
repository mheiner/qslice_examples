# Create Gibbs sampler to sample from the joint posterior
dhr <- function(n_iter, y, times, frequencies, period, m0, C0, W, a, b, burnin = 100, n_chains = 1, 
                sampler = "second_transform_pits", subtype = "",
                boundary = 0.01, width = 5, latent_scale = 10, parameters = NULL, verbose = TRUE, every = 1000) {
  
  require("coda")
  require("coda")
  require("truncnorm")
  require("crch")
  require("qslice")
  
  # Store the dhr function environment
  env <- environment()
  
  # Account for potential list of parameters
  if(!is.null(parameters)) {
    # If we have a list of parameters convert them to variables
    list2env(parameters, envir = env)
  }else{
    # If parameters were not supplied as a list, create a list to return
    parameters <- mget(ls(envir = env), envir = env)
  }
  
  # Create the source_all function
  source("scripts/ct_gibbs/source_all.R", local = TRUE)
  # Source all r scripts that apply to all samplers
  source_all("scripts/ct_gibbs/general", env)
  # Source all r scripts in the sampler folder, to the dhr environment
  source_all(paste0("scripts/ct_gibbs/", sampler), env)
  
  # Pull number of frequencies and length of series
  n <- length(y)
  n_freq <- length(frequencies)
  p <- n_freq + 1 # the plus 1 is for the drifting mean
  
  # Calculate the time difference to be used for scaling
  dift <- diff(times)
  time_diff <- c(mean(dift), dift)
  
  mcmc_list <- list()
  mcmc_list[["theta"]] <- vector("list", n_chains)
  mcmc_list[["sigma2"]] <- vector("list", n_chains)
  mcmc_list[["phases"]] <- vector("list", n_chains)

  # Prepare a vector to hold the rejections for each chain
  rejections <- numeric(n_chains)
  
  # Adjust posterior sample size to include burnin
  total_samples <- burnin + n_iter
  
  # Initialize objects to hold MCMC iterations
  theta_array <- array(dim = c(p, n, n_iter))
  sigma2 <- numeric(n_iter)
  phases <- matrix(nrow = n_freq, ncol = n_iter)

  # If we are using a PITS variant, capture the uniform variables
  pits_variant <- grepl("pits|Qslice", sampler)
  latent_variant <- grepl("latent", sampler)
  pgibbs_variant <- sampler_type == "pgibbs"
  if (isTRUE(pits_variant)) {
    # Storage for uniform transformed variables (for diagnostics on pseudo-targets)
    mcmc_list[["phi"]] <- vector("list", n_chains)
    phi_array <- array(dim = c(n_freq, n, n_iter))
  }

  for (chain_num in 1:n_chains) {
    
    ## Run a frequentist harmonic regression on corrupted data to get the 
    # starting phases and variance
    n_fourier <- 2 * n_freq # Allocate two fourier terms for each frequency
    regression_x <- matrix(1, nrow = n, ncol = n_fourier)
    # Populate the X matrix with the Fourier terms
    for(i in 1:n_fourier) {
      if((i %% 2) == 1){
        regression_x[,i] <- sin(2 * pi * frequencies[ceiling(i/2)] * times / period)
      } else {
        regression_x[,i] <- cos(2 * pi * frequencies[i/2] * times / period)
      }
    }
    # Run the regression
    y_corrupt <- y + rnorm(n, mean = 0, sd = 0.01*diff(range(y)))
    shr <- lm(y_corrupt ~ regression_x) # static harmonic regression with ols
    # Extract the fourier estimates to calculate phases
    fourier_estimates <- unname(coef(shr)[-1])
    # Obtain phase estimates
    sins <- fourier_estimates[seq(1, n_fourier, by = 2)] # Store the sin terms
    coss <- fourier_estimates[seq(2, n_fourier, by = 2)] # Store the cos terms
    
    # Extract the squared residual sd as a starting point for sigma2
    sigma2_now <- sigma(shr)^2
    
    phases_now <- atan2(-sins, coss) * period / (2 * pi) # Calculate phase and convert from (-pi,pi) to (0,period) for interpretability and to match the rest of the code
    phases_now <- phases_now + period * (phases_now < 0)
    
    # Establish known F values
    F_mat <- matrix(nrow = n, ncol = n_freq + 1)
    F_mat[,1] <- 1
    for(i in 1:n_freq) {
      F_mat[, i + 1] <- cos(2 * pi * (times * frequencies[i] + phases_now[i]) / period)
    }
    
    # Store initial value of theta from FFBS draw
    theta_now <- ffbs_rw(y = y, F_mat = F_mat, m0 = m0, C0 = C0, W = W, 
                         sigma2 = sigma2_now, time_diff = time_diff)
    
    theta_now[2:p,] <- pmax(theta_now[2:p,], 0.01) # Recognize that the intercept can go negative
    
    if (isTRUE(latent_variant)) {
      # Initialize s_mat for latent slice
      s_mat <- matrix(0.5*latent_scale, nrow = n_freq, ncol = n)
    } 
    # if (isTRUE(pgibbs_variant)) {
    #   theta0_now <- rtruncnorm(n_freq, mean = theta_now[2:p, 1], sd = sqrt(diag(C0)[2:p]/10), a = 0.0) # stopped sampling theta0 (init)
    # }
  
    # Initialize rejections
    chain_rejections <- 0
    
    # Begin MCMC iterations-----------------------------------------------------------
    
    for(j in 1:total_samples) {
      
      # Update user on sampling process as requested by the user
      if(verbose){
        if(j == burnin) {
          print(paste0("chain ", chain_num, " burnin of ", burnin, " complete"))
        } else if((j - burnin) %% every == 0 & (j - burnin) > 0) {
          print(paste0("chain ", chain_num, ", sample ", j - burnin, " of ", n_iter))
        }
      }
      
      ### Update all thetas
      tmp <- update_theta()
      theta_now <- tmp$theta
      chain_rejections <- chain_rejections + tmp$rejections
      if (isTRUE(pits_variant)) {
        phi_mat <- tmp$phi_mat
      }
      if (isTRUE(latent_variant)) {
        s_mat <- tmp$s_mat
      }
      # if (isTRUE(pgibbs_variant)) { # stopped sampling theta0 (init)
      #   theta0_now <- tmp$theta0
      # }

      ### Update the phase shifts
      phases_now <- phase_update(y, phases_now, sigma2_now, theta_now, frequencies, period)
      
      # Once the phase shifts are updated, update the F matrix
      for(i in 1:n_freq) {
        F_mat[, i + 1] <- cos(2 * pi * (times * frequencies[i] + phases_now[i]) / period)
      }
      
      ### Update sigma2
      sigma2_now <- sigma2_update(y, a, b, F_mat, theta_now)
      
      ### Fill in
      if (j > burnin) {
        jj <- j - burnin
        phases[,jj] <- phases_now
        sigma2[jj] <- sigma2_now
        theta_array[,,jj] <- theta_now
        if (isTRUE(pits_variant)) {
          phi_array[,,jj] <- phi_mat
        }
      }
      
      # cat(j, "\r")
    } # This is the end of the MCMC loop
    
    # Store the total rejections for that chain
    rejections[chain_num] <- chain_rejections
    
    mcmc_list[["theta"]][[chain_num]] <- theta_array
    mcmc_list[["sigma2"]][[chain_num]] <- sigma2
    mcmc_list[["phases"]][[chain_num]] <- phases
    if (isTRUE(pits_variant)) {
      mcmc_list[["phi"]][[chain_num]] <- phi_array
    }
  } # This is the end of the chains loop
  
  return_list <- list(
    mcmc_list = mcmc_list,
    parameters = parameters,
    rejections = rejections
  )
  
  return(
    return_list
  )
} # This is the end of the Gibbs function

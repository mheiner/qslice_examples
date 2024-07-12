
if (data_type == "river_long") {
  
  load("data/rivers_dhr.rda") # River data (see README)
  
  datc$river <- datc$StationC |> as.factor() |> as.numeric()
  river_data <- datc %>% filter(river == 1, time >= 1990, time < 2005) # this is the T = 117 series at Station 3052134
  
  stopifnot(p <= 4)
  C00 <- matrix(c(5^2, 0, 0, 0,
                  0, 4^2, 0, 0,
                  0, 0, 2^2, 0,
                  0, 0, 0, 2^2), nrow = 4, ncol = 4, byrow = TRUE)
  W0 <- matrix(12 * c(1/20^2, 0, 0, 0,
                      0, 1/8^2, 0, 0,
                      0, 0, 1/8^2, 0,
                      0, 0, 0, 1/8^2), nrow = 4, ncol = 4, byrow = TRUE)
  
  sample_parameters <- list(n_iter = n_iter,
                            y = river_data$lNO3[1:n],
                            frequencies = 1:(p-1),
                            period = 1,
                            m0 = c(0,0,0,0)[1:p],
                            C0 = C00[1:p, 1:p],
                            W = W0[1:p, 1:p],
                            a = 5/2,
                            b = 5*.25/2,
                            burnin = burnin,
                            n_chains = n_chains,
                            boundary = NA,
                            width = 1.5,
                            latent_scale = 0.3,
                            n_particles = 25,
                            pseu_family = "normal",
                            pseu_params = list(degf = Inf, sc_adj = 1.0),
                            times = river_data$time[1:n])
  
  data_sim <- river_data
  
} else if (data_type == "sim") {
  
  set.seed(1)
  
  tt <- seq(0, n/12, length = n)
  data_sim <- data.frame(time = tt + 1990,
                         lNO3 = 6 - tt + seq(0, 4.0, length = n) * sin(2*pi*tt) + rnorm(n, mean = 0, sd = 0.1))
  
  stopifnot(p <= 4)
  C00 <- matrix(c(5^2, 0, 0, 0,
                  0, 4^2, 0, 0,
                  0, 0, 2^2, 0,
                  0, 0, 0, 2^2), nrow = 4, ncol = 4, byrow = TRUE)
  W0 <- matrix(12 * c(1/20^2, 0, 0, 0,
                      0, 1/8^2, 0, 0,
                      0, 0, 1/8^2, 0,
                      0, 0, 0, 1/8^2), nrow = 4, ncol = 4, byrow = TRUE)
    
  sample_parameters <- list(n_iter = n_iter,
                            y = data_sim$lNO3[1:n],
                            frequencies = 1:(p-1),
                            period = 1,
                            m0 = c(0,0,0,0)[1:p],
                            C0 = C00[1:p, 1:p],
                            W = W0[1:p, 1:p],
                            a = 5/2,
                            b = 5*.25/2,
                            burnin = burnin,
                            n_chains = n_chains,
                            boundary = NA,
                            width = 1.5,
                            latent_scale = 0.3,
                            n_particles = 25,
                            pseu_family = "normal",
                            pseu_params = list(degf = Inf, sc_adj = 1.0),
                            times = data_sim$time[1:n])
}
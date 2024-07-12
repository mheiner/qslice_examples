
residual_analysis_DLM <- function(resids, fitted_values, bins = 15, confidence_level = 0.95) {
  library(ggplot2)
  library(patchwork)
  library(lomb)
  
  old_theme <- theme_set(theme_minimal())
  n <- length(resids)
  
  # Lomb-Scargle Periodogram
  perio <- lsp(resids, plot = FALSE, type = "period")
  perio_frame <- data.frame(period = perio$scanned, power = perio$power)
  label_frame <- data.frame(peak_at = perio$peak.at, peak = perio$peak, label = round(perio$peak.at[1],1))
  p1 <- ggplot(perio_frame, aes(period, power))
  p1 <- p1 + geom_line() + geom_hline(yintercept = perio$sig.level,
                                      lty = 2,
                                      color = "blue")
  p1 <- p1 + labs(x = "Period", y = "Power")
  p1 <- p1 + geom_text(mapping = aes(x = peak_at[1], y = peak, label = label),
                       data = label_frame,
                       hjust = -0.5,
                       vjust = 0,
                       inherit.aes = FALSE)

  # Fitted values vs. residuals
  p2 <- ggplot(mapping = aes(fitted_values, resids))
  p2 <- p2 + geom_point()
  p2 <- p2 + labs(x = "Fitted Values", y = "Residuals")
  
  # Residuals over time
  p3 <- ggplot(mapping = aes(seq_len(n), resids))
  p3 <- p3 + labs(x = "Time", y = "Residuals")
  p3 <- p3 + geom_line()
  
  # ACF data
  # resids.acf <- acf(resids, lag.max = 24, plot = FALSE)
  # ## ggplot() of ACF function
  # acf.frame <- data.frame(Lag = resids.acf$lag, ACF = resids.acf$acf)
  # ggplot(data = acf.frame, aes(x = Lag, y = ACF)) + geom_col()
  
  # PACF
  resids_pacf <- pacf(resids, lag.max = 24, plot = FALSE)
  ## ggplot() of PACF function
  pacf_frame <- data.frame(Lag = resids_pacf$lag, PACF = resids_pacf$acf)
  interval <- c(
    upper = qnorm((1 + confidence_level)/2)/sqrt(resids_pacf$n.used),
    lower = -qnorm((1 + confidence_level)/2)/sqrt(resids_pacf$n.used)
  )
  
  p4 <- ggplot(data = pacf_frame, aes(x = Lag, y = PACF)) + geom_col()
  p4 <- p4 + geom_hline(yintercept = interval, col = "blue", lty = 2)
  p4 <- p4 + labs(y = "Partial ACF")
  
  resids_frame <- data.frame("Residuals" = resids)
  density_frame <- data.frame("x" = seq(min(resids), max(resids), length = 1000),
                              "y" = dnorm(seq(min(resids), max(resids), length = 1000), mean(resids), sd(resids)))
  p5 <- ggplot(data = resids_frame, aes(x = Residuals)) + geom_histogram(aes(y = after_stat(density)), bins = bins) + geom_line(data = density_frame, mapping = aes(x, y), col = "red", linewidth = 1)
  p5 <- p5 + labs(y = "Density")
  
  theme_set(old_theme)
  
  print(
    (p1 + p2)/(p5 + p4)
  )
  
}

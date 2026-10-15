## suggested by ChatGPT Codex to replace package "crch"

.ptt0 <- function(q, location, scale, df) {
  log_surv_0 <- pt(
    (0 - location) / scale,
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )

  log_surv_q <- pt(
    (q - location) / scale,
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )

  # 1 - S(q) / S(0), evaluated accurately when q is near zero
  out <- -expm1(log_surv_q - log_surv_0)

  below <- !is.na(q) & q < 0
  out[below] <- 0

  out
}

.qtt0 <- function(p, location, scale, df) {
  log_surv_0 <- pt(
    (0 - location) / scale,
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )

  log_surv_q <- log1p(-p) + log_surv_0

  location +
    scale *
      qt(
        log_surv_q,
        df = df,
        lower.tail = FALSE,
        log.p = TRUE
      )
}

.ptnorm0 <- function(q, location, scale) {
  log_surv_0 <- pnorm(
    (0 - location) / scale,
    lower.tail = FALSE,
    log.p = TRUE
  )

  log_surv_q <- pnorm(
    (q - location) / scale,
    lower.tail = FALSE,
    log.p = TRUE
  )

  # 1 - S(q) / S(0), evaluated accurately when q is near zero
  out <- -expm1(log_surv_q - log_surv_0)

  below <- !is.na(q) & q < 0
  out[below] <- 0

  out
}

.qtnorm0 <- function(p, location, scale) {
  log_surv_0 <- pnorm(
    (0 - location) / scale,
    lower.tail = FALSE,
    log.p = TRUE
  )

  log_surv_q <- log1p(-p) + log_surv_0

  location +
    scale *
      qnorm(
        log_surv_q,
        lower.tail = FALSE,
        log.p = TRUE
      )
}

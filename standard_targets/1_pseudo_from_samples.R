if (target == "igamma") {
  df_samples <- 1
} else if (target == "igammalog") {
  df_samples <- c(1, 5)
} else if (target == "gamma") {
  df_samples <- c(1, 5)
} else {
  df_samples <- c(1, 5, 20)
}

init_fn <- function() {
  c(
    loc = rnorm(1, mean = median(samples_train), sd = 0.1),
    sc = abs(rnorm(
      1,
      mean = diff(quantile(samples_train, probs = c(0.25, 0.75))),
      sd = 0.1
    ))
  )
}

ps_success <- FALSE
ps_attempts <- 0L
max_ps_attempts <- 5L

if (isTRUE(run_info$subtype == "MSW_samples")) {
  while (!ps_success && (ps_attempts < max_ps_attempts)) {
    ps_attempts <- ps_attempts + 1L
    new_pseu <- tryCatch(
      {
        out <- pseudo_opt(
          samples = samples_train,
          type = "samples",
          family = "t",
          df = df_samples,
          lb = truth$lb,
          ub = truth$ub,
          init_fn = init_fn,
          utility_type = "MSW",
          plot = FALSE
        )
        ps_success <- TRUE
        out
      },
      error = function(e) {
        if (
          grepl("All pseudo-target utility optimizations failed.", e$message)
        ) {
          message(sprintf(
            "Attempt %d failed with message: %s. Retrying...",
            ps_attempts,
            e$message
          ))
          NULL
        } else {
          stop(e)
        }
      }
    )
  }
  if (!ps_success) {
    stop("All rounds of pseudo-target optimization failed.")
  }

  trials[["Qslice"]][["MSW_samples"]][["pseudo"]] <- new_pseu
} else if (isTRUE(run_info$subtype == "AUC_samples")) {
  while (!ps_success && (ps_attempts < max_ps_attempts)) {
    ps_attempts <- ps_attempts + 1L
    new_pseu <- tryCatch(
      {
        out <- pseudo_opt(
          samples = samples_train,
          type = "samples",
          family = "t",
          df = df_samples,
          lb = truth$lb,
          ub = truth$ub,
          init_fn = init_fn,
          utility_type = "AUC",
          plot = FALSE
        )
        ps_success <- TRUE
        out
      },
      error = function(e) {
        if (
          grepl("All pseudo-target utility optimizations failed.", e$message)
        ) {
          message(sprintf(
            "Attempt %d failed with message: %s. Retrying...",
            ps_attempts,
            e$message
          ))
          NULL
        } else {
          stop(e)
        }
      }
    )
  }
  if (!ps_success) {
    stop("All rounds of pseudo-target optimization failed.")
  }

  trials[["Qslice"]][["AUC_samples"]][["pseudo"]] <- new_pseu
} else if (isTRUE(run_info$subtype == "MM_Cauchy")) {
  new_pseu <- list(
    # for uniformity of structure
    pseu = pseudo_list(
      family = "cauchy",
      params = list(loc = mean(samples_train), sc = sd(samples_train)),
      lb = truth$lb,
      ub = truth$ub
    )
  )

  trials[["Qslice"]][["MM_Cauchy"]][["pseudo"]] <- new_pseu
}

trials[[run_info$type]][[run_info$subtype]]$algo_descrip <- trials[[
  run_info$type
]][[run_info$subtype]]$pseudo$pseu$txt

rm(new_pseu)

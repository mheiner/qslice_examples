Qtypes <- c(
  "MSW",
  "MSW_samples",
  "AUC",
  "AUC_wide",
  "AUC_samples",
  "MM_Cauchy",
  "Laplace_Cauchy"
)

trials[["Qslice"]] <- list()

trials[["Qslice"]][["MSW"]] <- list(
  target = target,
  type = "Qslice",
  subtype = "MSW",
  pseudo = pseudo_opt(
    log_target = truth$ld,
    type = "function",
    family = "t",
    lb = truth$lb,
    ub = truth$ub,
    utility_type = "MSW",
    plot = TRUE
  )
)

trials[["Qslice"]][["MSW_samples"]] <- list(
  target = target,
  type = "Qslice",
  subtype = "MSW_samples",
  pseudo = list(pseu = list(txt = "TBD"))
)

trials[["Qslice"]][["AUC"]] <- list(
  target = target,
  type = "Qslice",
  subtype = "AUC",
  pseudo = pseudo_opt(
    log_target = truth$ld,
    type = "function",
    family = "t",
    lb = truth$lb,
    ub = truth$ub,
    utility_type = "AUC",
    plot = TRUE
  )
)

trials[["Qslice"]][["AUC_wide"]] <- list(
  target = target,
  type = "Qslice",
  subtype = "AUC_wide",
  pseudo = list(
    # for uniformity of structure
    pseu = pseudo_list(
      family = "t",
      params = list(
        loc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$loc,
        sc = wide_factor * trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$sc,
        degf = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$degf
      ),
      lb = truth$lb,
      ub = truth$ub
    )
  )
)

trials[["Qslice"]][["AUC_samples"]] <- list(
  target = target,
  type = "Qslice",
  subtype = "AUC_samples",
  pseudo = list(pseu = list(txt = "TBD"))
)

trials[["Qslice"]][["MM_Cauchy"]] <- list(
  target = target,
  type = "Qslice",
  subtype = "MM_Cauchy",
  pseudo = list(pseu = list(txt = "TBD"))
)

trials[["Qslice"]][["Laplace_Cauchy"]] <- list(
  target = target,
  type = "Qslice",
  subtype = "Laplace_Cauchy",
  pseudo = list(
    # for uniformity of structure
    pseu = lapprox(
      log_target = truth$d,
      init = init_lapprox,
      family = "cauchy",
      sc_adj = 1.0,
      lb = truth$lb,
      ub = truth$ub
    )
  )
)

stopifnot(all.equal(Qtypes, names(trials[["Qslice"]])))
for (xx in Qtypes) {
  trials[["Qslice"]][[xx]][["algo_descrip"]] <- trials[["Qslice"]][[
    xx
  ]]$pseudo$pseu$txt
}


## Independence M-H

trials[["imh"]] <- list()

trials[["imh"]][["AUC"]] <- list(
  target = target,
  type = "imh",
  subtype = "AUC",
  pseudo = list(
    # for uniformity of structure
    pseu = pseudo_list(
      family = "t",
      params = list(
        loc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$loc,
        sc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$sc,
        degf = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$degf
      ),
      lb = truth$lb,
      ub = truth$ub
    )
  )
)

trials[["imh"]][["AUC_wide"]] <- list(
  target = target,
  type = "imh",
  subtype = "AUC_wide",
  pseudo = list(
    # for uniformity of structure
    pseu = pseudo_list(
      family = "t",
      params = list(
        loc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$loc,
        sc = wide_factor * trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$sc,
        degf = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$degf
      ),
      lb = truth$lb,
      ub = truth$ub
    )
  )
)

for (xx in names(trials[["imh"]])) {
  trials[["imh"]][[xx]][["algo_descrip"]] <- trials[["imh"]][[
    xx
  ]]$pseudo$pseu$txt
}


## Generalized elliptical slice

trials[["gess"]] <- list()

trials[["gess"]][["AUC"]] <- list(
  target = target,
  type = "gess",
  subtype = "AUC",
  pseudo = list(
    loc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$loc,
    sc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$sc,
    degf = trials[["Qslice"]][["AUC"]]$pseudo$pseu$params$degf
  )
)

for (xx in names(trials[["gess"]])) {
  trials[["gess"]][[xx]][["algo_descrip"]] <- trials[["Qslice"]][[
    xx
  ]]$pseudo$pseu$txt
}


## create schedule
sched0[["Qslice"]] <- expand.grid(
  target = target,
  type = "Qslice",
  subtype = Qtypes,
  stringsAsFactors = FALSE
)
sched0[["Qslice"]]$algo_descrip <- sapply(
  sched0[["Qslice"]]$subtype,
  function(xx) {
    trials[["Qslice"]][[xx]]$algo_descrip
  }
)

sched0[["Qslice"]]

sched0[["imh"]] <- expand.grid(
  target = target,
  type = "imh",
  subtype = names(trials[["imh"]]),
  stringsAsFactors = FALSE
)
sched0[["imh"]]$algo_descrip <- sapply(sched0[["imh"]]$subtype, function(xx) {
  trials[["imh"]][[xx]]$algo_descrip
})

sched0[["imh"]]


sched0[["gess"]] <- expand.grid(
  target = target,
  type = "gess",
  subtype = names(trials[["gess"]]),
  stringsAsFactors = FALSE
)
sched0[["gess"]]$algo_descrip <- sapply(sched0[["gess"]]$subtype, function(xx) {
  trials[["gess"]][[xx]]$algo_descrip
})

sched0[["gess"]]

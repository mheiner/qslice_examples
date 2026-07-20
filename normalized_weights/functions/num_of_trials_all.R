# how many lines are in a rds file
# author: Sam Johnson

args <- commandArgs(trailingOnly = TRUE)

dte <- args[1] |> as.numeric()

file <- paste0("input/schedule_all_", dte, ".rda")

load(file)

cat(n_jobs)

quit(save = "no")

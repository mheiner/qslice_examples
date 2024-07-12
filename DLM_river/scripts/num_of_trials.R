# how many lines are in a rds file
# authors: Sam Johnson, Matt Heiner

args <- commandArgs(trailingOnly = TRUE)
data_type <- args[1]
round <- args[2]

file <- paste0("data/schedule_data", data_type, "_", round, ".rda")

load(file)

cat(n_jobs)

quit(save = "no")

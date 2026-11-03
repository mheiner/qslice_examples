dte <- 261019

files_all <- list.files("output")
files_use <- files_all[grep(paste0(".*_dte", dte), x = files_all)]

temp <- lapply(paste0("output/", files_use), read.csv)
out <- do.call(rbind, temp)
str(out)

## clean up
# system(paste("rm ", paste0("output/*.csv")))
# system(paste("rm ", paste0("logs/*.txt")))

(file_write <- paste0("output/combined_all_", dte, ".csv"))
write.csv(
  file = file_write,
  out,
  row.names = FALSE
)

# system("ls output")

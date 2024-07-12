source_all <- function(folder, env) {
  r_scripts <- list.files(path = folder, pattern = "*.R")
  for(script in r_scripts) {
    script_path <- paste0(folder, "/", script)
    source(script_path, local = env)
  }
}

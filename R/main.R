# Entry point: Rscript /app/R/main.R
# Reads INPUT_DIR / OUTPUT_DIR and the parameters from the environment (see app.yml).

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  app_r <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else Sys.getenv("APP_R_DIR", "R")
  for (f in c("common.R", "fcs_io.R", "cytonorm.R")) sys.source(file.path(app_r, f), envir = globalenv())
  assign("APP_R_DIR", app_r, envir = globalenv())
})

main <- function() {
  input_dir <- Sys.getenv("INPUT_DIR", "/data/input")
  output_dir <- Sys.getenv("OUTPUT_DIR", "/data/output")
  status <- tryCatch({
    p <- cytonorm_params()
    run_cytonorm(input_dir, output_dir, APP_R_DIR, p)
    0L
  }, error = function(e) {
    log_line("ERROR: ", conditionMessage(e))
    1L
  })
  quit(save = "no", status = status)
}

if (sys.nframe() == 0L) main()

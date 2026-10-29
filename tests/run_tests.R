# Run the unit tests: Rscript tests/run_tests.R (from the repository root).
# Needs R with the packages the image installs (flowCore, FlowSOM, CytoNorm, ggplot2,
# jsonlite); `make test` runs this in the image when no local R is available.

suppressPackageStartupMessages({
  library(flowCore)
  library(FlowSOM)
  library(CytoNorm)
})
for (f in c("R/common.R", "R/fcs_io.R", "R/cytonorm.R", "tests/helpers.R", "tests/test_cytonorm.R")) source(f)
sink_log <- function(code) {
  # keep the per-sample log lines out of the test output
  out <- utils::capture.output(res <- tryCatch({ force(code); NULL }, error = function(e) e), type = "output")
  if (!is.null(res)) {
    cat(utils::tail(out, 15), sep = "\n")
    stop(res)
  }
  invisible(res)
}

tests <- ls(pattern = "^test_")
failed <- 0
for (t in tests) {
  t0 <- Sys.time()
  ok <- tryCatch({ sink_log(suppressWarnings(suppressMessages(get(t)()))); TRUE },
                 error = function(e) { cat("FAIL ", t, ": ", conditionMessage(e), "\n", sep = ""); FALSE })
  if (ok) cat(sprintf("ok   %s (%.1fs)\n", t, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  failed <- failed + !ok
}
cat(sprintf("%d test(s), %d failed\n", length(tests), failed))
quit(status = if (failed) 1 else 0)

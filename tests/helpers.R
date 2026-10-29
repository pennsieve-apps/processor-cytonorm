# Test helpers: a tiny synthetic fixture that follows the file contract, and assertions.

expect_true <- function(x, msg = deparse(substitute(x))) {
  if (!isTRUE(x)) stop("expectation failed: ", msg, call. = FALSE)
}
expect_equal_num <- function(a, b, tol = 1e-6, msg = "") {
  d <- max(abs(as.numeric(a) - as.numeric(b)))
  if (!is.finite(d) || d > tol) stop(sprintf("not equal (max diff %g > %g) %s", d, tol, msg), call. = FALSE)
}
expect_error <- function(expr, pattern = NULL) {
  err <- tryCatch({ force(expr); NULL }, error = function(e) e)
  if (is.null(err)) stop("expected an error", call. = FALSE)
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(err))) {
    stop("error did not match '", pattern, "': ", conditionMessage(err), call. = FALSE)
  }
  invisible(err)
}

FIX_CHANNELS <- c("FSC-A", "SSC-A", "BV421-A", "PE-A", "APC-A", "FITC-A", "AF-A", "Time")
FIX_MARKERS <- c("", "", "CD3", "CD19", "CD45RA", "GZMB", "", "")

#' Synthetic samples: two cell populations (CD3+ T cells / CD19+ B cells), marker values on an
#' arcsinh scale, and a per-batch technical effect (marker shift and stretch) that the same
#' "control" sample shows in every batch. Scatter, AF and Time are left untouched by
#' normalisation and are used to check that.
make_synthetic <- function(dir, batches = c("b1", "b2", "b3"), n = 3000, seed = 1,
                           extra_sample = NULL, missing_ref_batch = NULL) {
  set.seed(seed)
  dir.create(file.path(dir, "fcs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir, "metadata"), recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  for (bi in seq_along(batches)) {
    b <- batches[bi]
    shift <- c(0, 0.25, -0.2, 0.15, 0.3, -0.1)[bi]
    stretch <- c(1, 1.1, 0.92, 1.05, 0.95, 1.08)[bi]
    for (s in c("ctrl", "s1", "s2")) {
      if (identical(b, missing_ref_batch) && s == "ctrl") next
      t_frac <- switch(s, ctrl = 0.6, s1 = 0.7, s2 = 0.5)
      cond <- if (s == "s2") "stim" else "unstim"
      nt <- round(n * t_frac)
      pop_t <- c(rep(TRUE, nt), rep(FALSE, n - nt))
      cd3 <- ifelse(pop_t, rnorm(n, 3, 0.35), rnorm(n, 0.3, 0.3))
      cd19 <- ifelse(pop_t, rnorm(n, 0.2, 0.3), rnorm(n, 2.8, 0.4))
      cd45ra <- ifelse(pop_t, rnorm(n, 1.5, 0.8), rnorm(n, 2.5, 0.4))
      gzmb <- ifelse(pop_t, rnorm(n, 1 + (s == "s2"), 0.5), rnorm(n, 0.3, 0.3))
      mk <- cbind(cd3, cd19, cd45ra, gzmb) * stretch + shift
      m <- cbind(rnorm(n, 1e5, 2e4), rnorm(n, 5e4, 1e4), mk, rnorm(n, 300, 50), sort(runif(n, 0, 600)))
      colnames(m) <- FIX_CHANNELS
      sid <- paste0(b, "_", s)
      kw <- c("$CYT" = "Synthetic", "$DATE" = "01-Jan-2024", "$FIL" = paste0(sid, ".fcs"),
              "GROUPNAME" = "donor A", "PENNSIEVE_PREP" = "transform=arcsinh;cofactor=6000")
      for (i in seq_along(FIX_CHANNELS)) {
        kw[sprintf("$P%dN", i)] <- FIX_CHANNELS[i]
        kw[sprintf("$P%dR", i)] <- "262144"
        if (nzchar(FIX_MARKERS[i])) kw[sprintf("$P%dS", i)] <- FIX_MARKERS[i]
      }
      write_fcs(file.path(dir, "fcs", paste0(sid, ".fcs")), m, kw)
      rows[[length(rows) + 1]] <- data.frame(sample_id = sid, file = paste0(sid, ".fcs"), batch = b,
                                             condition = cond, n_events = n, control = s,
                                             stringsAsFactors = FALSE)
    }
  }
  samples <- do.call(rbind, rows)
  if (!is.null(extra_sample)) {
    file.copy(file.path(dir, "fcs", paste0(samples$sample_id[1], ".fcs")),
              file.path(dir, "fcs", paste0(extra_sample, ".fcs")))
  }
  write.csv(samples, file.path(dir, "metadata", "samples.csv"), row.names = FALSE)
  type <- ifelse(nzchar(FIX_MARKERS), "marker", ifelse(FIX_CHANNELS == "Time", "time",
                 ifelse(FIX_CHANNELS == "AF-A", "other", "scatter")))
  panel <- data.frame(channel = FIX_CHANNELS, marker = ifelse(nzchar(FIX_MARKERS), FIX_MARKERS, FIX_CHANNELS),
                      type = type, cofactor = ifelse(type == "marker", "6000", ""),
                      clustering = ifelse(FIX_MARKERS %in% c("CD3", "CD19", "CD45RA"), "true", "false"),
                      state = ifelse(FIX_MARKERS == "GZMB", "true", "false"))
  write.csv(panel, file.path(dir, "metadata", "panel.csv"), row.names = FALSE)
  invisible(samples)
}

with_env <- function(vars, code) {
  old <- Sys.getenv(names(vars), unset = NA, names = TRUE)
  do.call(Sys.setenv, as.list(vars))
  on.exit({
    for (k in names(old)) if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, stats::setNames(list(old[[k]]), k))
  })
  force(code)
}

new_tmp <- function(name) {
  d <- file.path(tempdir(), paste0(name, "-", as.integer(runif(1, 1, 1e9))))
  dir.create(d, recursive = TRUE)
  d
}

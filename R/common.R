# Shared helpers for the immune cytometry R apps: parameters, input discovery, metadata,
# logging and provenance. See workflow-immune-cytometry/docs/conventions.md.

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(...), "\n", sep = "")
  flush(stdout())
}

# ---- parameters (environment variables named as in app.yml) -------------------------

env_str <- function(name, default = "") {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v) || !nzchar(trimws(v))) default else trimws(v)
}

env_int <- function(name, default, min = -Inf, max = Inf) {
  v <- env_str(name, NA)
  if (is.na(v)) return(as.integer(default))
  x <- suppressWarnings(as.numeric(v))
  if (is.na(x) || x != round(x)) stop(sprintf("%s must be an integer, got '%s'", name, v))
  if (x < min || x > max) stop(sprintf("%s must be between %s and %s, got %s", name, min, max, v))
  as.integer(x)
}

env_choice <- function(name, default, choices) {
  v <- env_str(name, default)
  if (!tolower(v) %in% tolower(choices)) {
    stop(sprintf("%s must be one of %s, got '%s'", name, paste(choices, collapse = ", "), v))
  }
  choices[match(tolower(v), tolower(choices))]
}

env_list <- function(name) {
  v <- env_str(name, "")
  if (!nzchar(v)) return(character())
  x <- trimws(strsplit(v, ",", fixed = TRUE)[[1]])
  x[nzchar(x)]
}

# ---- inputs ---------------------------------------------------------------------------

# All files under `dir` (symlinks followed), as absolute paths.
list_all_files <- function(dir) {
  if (!dir.exists(dir)) return(character())
  f <- list.files(dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  f[!dir.exists(f)]
}

#' Find the FCS files: INPUT_DIR/<port>/ when that folder exists (linked input), else
#' every *.fcs anywhere under INPUT_DIR (unlinked run: everything merged at the root).
find_fcs <- function(input_dir, port = "fcs") {
  sub <- file.path(input_dir, port)
  files <- if (dir.exists(sub)) list_all_files(sub) else list_all_files(input_dir)
  files <- files[grepl("\\.fcs$", files, ignore.case = TRUE)]
  sort(files)
}

#' Find a metadata file by name: INPUT_DIR/<port>/<name> first, else the shallowest file
#' with that name anywhere under INPUT_DIR.
find_named <- function(input_dir, name, port = "metadata") {
  direct <- file.path(input_dir, port, name)
  if (file.exists(direct)) return(direct)
  files <- list_all_files(input_dir)
  hits <- files[basename(files) == name]
  if (!length(hits)) return(NA_character_)
  depth <- lengths(strsplit(hits, "/", fixed = TRUE))
  hits[order(depth, hits)][1]
}

rel_path <- function(path, root) {
  root <- sub("/+$", "", normalizePath(root, mustWork = FALSE))
  p <- path
  ifelse(startsWith(p, paste0(root, "/")), substring(p, nchar(root) + 2), p)
}

read_csv_chr <- function(path) {
  utils::read.csv(path, colClasses = "character", check.names = FALSE, na.strings = character(),
                  stringsAsFactors = FALSE, strip.white = TRUE, fileEncoding = "UTF-8")
}

#' samples.csv: one row per sample; sample_id required, plus whatever columns `needed`.
read_samples <- function(path, needed = character()) {
  s <- read_csv_chr(path)
  miss <- setdiff(c("sample_id", needed), names(s))
  if (length(miss)) stop(sprintf("%s lacks column(s): %s", path, paste(miss, collapse = ", ")))
  if (anyDuplicated(s$sample_id)) stop(sprintf("%s has duplicate sample_id values", path))
  s
}

#' panel.csv: channel, marker, type, cofactor, clustering, state.
read_panel <- function(path) {
  p <- read_csv_chr(path)
  miss <- setdiff(c("channel", "marker", "type"), names(p))
  if (length(miss)) stop(sprintf("%s lacks column(s): %s", path, paste(miss, collapse = ", ")))
  for (col in c("clustering", "state")) {
    p[[col]] <- if (col %in% names(p)) tolower(p[[col]]) %in% c("true", "1", "yes") else FALSE
  }
  if (!"cofactor" %in% names(p)) p$cofactor <- ""
  p
}

#' The marker channels to work on: MARKERS (comma-separated marker or channel names) or,
#' when empty, every panel row of type marker. Returns the panel rows, in panel order.
select_markers <- function(panel, requested = character()) {
  if (!length(requested)) {
    sel <- panel[panel$type == "marker", , drop = FALSE]
  } else {
    hit <- panel$marker %in% requested | panel$channel %in% requested
    unknown <- requested[!requested %in% c(panel$marker, panel$channel)]
    if (length(unknown)) stop(sprintf("MARKERS not in panel.csv: %s", paste(unknown, collapse = ", ")))
    sel <- panel[hit, , drop = FALSE]
  }
  if (!nrow(sel)) stop("no marker channels selected (panel.csv has no rows of type 'marker'?)")
  sel
}

#' Match FCS files to samples.csv rows by sample_id (= file name without .fcs).
match_samples <- function(fcs_files, samples) {
  ids <- sub("\\.fcs$", "", basename(fcs_files), ignore.case = TRUE)
  dup <- ids[duplicated(ids)]
  if (length(dup)) stop(sprintf("two input FCS files share a sample_id: %s", paste(unique(dup), collapse = ", ")))
  data.frame(sample_id = ids, path = fcs_files,
             row = match(ids, samples$sample_id), stringsAsFactors = FALSE)
}

# ---- provenance -------------------------------------------------------------------------

app_version <- function(app_dir) {
  for (f in c(file.path(app_dir, "VERSION"), file.path(dirname(app_dir), "VERSION"), "/app/VERSION")) {
    if (file.exists(f)) return(trimws(readLines(f, n = 1, warn = FALSE)))
  }
  "unknown"
}

package_versions <- function(pkgs) {
  v <- lapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA))
  names(v) <- pkgs
  c(list(R = paste(R.version$major, R.version$minor, sep = ".")), v)
}

write_json <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
}

write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

iso_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

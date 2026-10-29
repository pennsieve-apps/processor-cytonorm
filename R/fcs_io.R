# FCS reading and writing that keeps a file's keywords.
#
# flowCore::write.FCS rewrites a file's TEXT segment (adds flowCore_$PnRmin/max, changes
# $PnR, writes FCS3.0), so files are written here directly: FCS 3.1, float32
# little-endian, list mode, with every keyword of the input file kept except the
# segment offsets and the per-channel data-format keywords that a float32 file must
# have ($PnB = 32, $PnE = 0,0). Data are read with flowCore (no transformation, no
# range truncation), keywords with read.FCSheader (the TEXT segment as written).

SEGMENT_KEYWORDS <- c("$BEGINANALYSIS", "$ENDANALYSIS", "$BEGINSTEXT", "$ENDSTEXT",
                      "$BEGINDATA", "$ENDDATA", "$NEXTDATA")

# Keywords flowCore adds when reading that are not in the file.
FLOWCORE_ADDED <- c("FCSversion", "FILENAME", "GUID", "ORIGINALGUID", "transformation",
                    "$CYTOLIB_VERSION", "CREATOR_CYTOLIB")

#' Read one FCS file: list(exprs = numeric matrix (events x channels, colnames $PnN),
#' keywords = named character vector of the TEXT segment, markers = $PnS per channel).
read_fcs <- function(path, which_lines = NULL) {
  kw <- flowCore::read.FCSheader(path, emptyValue = FALSE)[[1]]
  kw <- kw[!names(kw) %in% FLOWCORE_ADDED]
  ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE,
                           emptyValue = FALSE, which.lines = which_lines)
  m <- flowCore::exprs(ff)
  pd <- flowCore::pData(flowCore::parameters(ff))
  markers <- as.character(pd$desc)
  markers[is.na(markers)] <- ""
  names(markers) <- as.character(pd$name)
  list(exprs = m, keywords = kw, markers = markers, flowframe = ff)
}

#' Keyword lookup, case-insensitive (FCS keywords are case-insensitive).
fcs_keyword <- function(keywords, key, default = NA_character_) {
  i <- match(toupper(key), toupper(names(keywords)))
  if (is.na(i)) default else unname(keywords[[i]])
}

.fcs_text <- function(kw, delim) {
  esc <- function(x) gsub(delim, paste0(delim, delim), x, fixed = TRUE)
  body <- paste0(esc(names(kw)), delim, esc(unname(kw)), delim, collapse = "")
  paste0(delim, body)
}

#' Write a matrix as an FCS 3.1 float32 file.
#'
#' @param path output path
#' @param exprs numeric matrix, events x channels; colnames must be the $PnN of `keywords`
#'   in the same order (the channel set of the source file)
#' @param keywords named character vector of the source file's keywords
#' @param extra named character vector of keywords to add or replace (e.g. PENNSIEVE_NORM)
write_fcs <- function(path, exprs, keywords, extra = character()) {
  stopifnot(is.matrix(exprs))
  n <- nrow(exprs)
  p <- ncol(exprs)
  kw <- as.list(keywords)
  up <- toupper(names(kw))
  # Channel names in the keywords must match the matrix columns.
  pn <- vapply(seq_len(p), function(i) {
    v <- kw[[match(sprintf("$P%dN", i), up)]]
    if (is.null(v)) NA_character_ else v
  }, character(1))
  if (anyNA(pn) || !identical(unname(pn), unname(as.character(colnames(exprs))))) {
    stop("write_fcs: matrix columns do not match the $PnN keywords of the source file")
  }
  drop <- up %in% SEGMENT_KEYWORDS | up %in% c("$TOT", "$PAR", "$DATATYPE", "$BYTEORD", "$MODE") |
    grepl("^\\$P[0-9]+[BE]$", up)
  kw <- kw[!drop]
  for (k in names(extra)) {
    j <- match(toupper(k), toupper(names(kw)))
    if (!is.na(j)) kw[[j]] <- extra[[k]] else kw[[k]] <- extra[[k]]
  }
  fixed <- c("$BYTEORD" = "1,2,3,4", "$DATATYPE" = "F", "$MODE" = "L",
             "$PAR" = as.character(p), "$TOT" = as.character(n),
             "$NEXTDATA" = "0", "$BEGINANALYSIS" = "0", "$ENDANALYSIS" = "0",
             "$BEGINSTEXT" = "0", "$ENDSTEXT" = "0")
  for (i in seq_len(p)) {
    fixed[sprintf("$P%dB", i)] <- "32"
    fixed[sprintf("$P%dE", i)] <- "0,0"
  }
  kw <- c(as.list(fixed), kw)
  vals <- vapply(kw, function(v) as.character(v)[1], character(1))
  vals <- trimws(vals)
  vals <- vals[!is.na(vals) & nzchar(vals)]  # FCS 3.1 does not allow empty values
  delim <- "|"
  data_bytes <- as.numeric(n) * p * 4
  text_start <- 58
  begin <- 0
  repeat {
    end <- if (data_bytes > 0) begin + data_bytes - 1 else begin
    all_kw <- c(vals, "$BEGINDATA" = sprintf("%.0f", begin), "$ENDDATA" = sprintf("%.0f", end))
    text <- enc2utf8(.fcs_text(all_kw, delim))
    text_bytes <- nchar(text, type = "bytes")
    new_begin <- text_start + text_bytes
    if (new_begin == begin) break
    begin <- new_begin
  }
  hdr_off <- function(x) if (x > 99999999) sprintf("%8s", "0") else sprintf("%8.0f", x)
  header <- paste0("FCS3.1", "    ",
                   sprintf("%8.0f", text_start), sprintf("%8.0f", text_start + text_bytes - 1),
                   hdr_off(begin), hdr_off(end), sprintf("%8s", "0"), sprintf("%8s", "0"))
  stopifnot(nchar(header, type = "bytes") == 58)
  con <- file(path, "wb")
  on.exit(close(con))
  writeBin(charToRaw(header), con)
  writeBin(charToRaw(text), con)
  chunk <- 200000L
  if (n > 0) {
    for (s in seq(1L, n, by = chunk)) {
      e <- min(n, s + chunk - 1L)
      writeBin(as.vector(t(exprs[s:e, , drop = FALSE])), con, size = 4, endian = "little")
    }
  }
  invisible(path)
}

# processor-cytonorm: CytoNorm batch normalisation (Van Gassen et al. 2020; CytoNorm 2.0,
# Quintelier et al. 2025) of arcsinh-transformed FCS files.
#
# 1. Training data, one aggregate per batch:
#    MODE=controls  the reference samples of the batch (REFERENCE_COLUMN in REFERENCE_VALUES,
#                   e.g. the same control aliquot measured in every batch) - CytoNorm 1.
#    MODE=cohort    every sample of the batch - CytoNorm 2.0's "aggregate as a proxy for a
#                   control" (vignette "CytoNorm without controls"); assumes the batches hold
#                   comparable samples (same mix of conditions).
#    At most TRAIN_CELLS_PER_SAMPLE events per sample (seeded), so memory stays bounded.
# 2. CytoNorm.train: FlowSOM on the panel's clustering markers (N_METACLUSTERS metaclusters),
#    then per metacluster, batch and marker N_QUANTILES quantiles and a monotone spline to
#    the goal quantiles (GOAL = mean over batches, or one batch).
# 3. Every sample (references included) is mapped to the FlowSOM model and the spline of its
#    batch and metacluster applied to each marker - the same arithmetic as
#    CytoNorm.normalize, done in memory one file at a time, so a file's other channels and
#    keywords are written back unchanged. Data are already arcsinh-transformed, so the
#    transform lists are identity (NULL).

CYTONORM_PKGS <- c("CytoNorm", "FlowSOM", "flowCore", "ggplot2", "jsonlite")
QC_PROBS <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)

cytonorm_params <- function() {
  mode <- env_choice("MODE", "controls", c("controls", "cohort"))
  p <- list(
    MODE = mode,
    BATCH_COLUMN = env_str("BATCH_COLUMN", "batch"),
    REFERENCE_COLUMN = env_str("REFERENCE_COLUMN", "control"),
    REFERENCE_VALUES = env_list("REFERENCE_VALUES"),
    GOAL = env_str("GOAL", "mean"),
    N_METACLUSTERS = env_int("N_METACLUSTERS", 5, min = 3, max = 100),
    N_QUANTILES = env_int("N_QUANTILES", 101, min = 2, max = 1001),
    XDIM = env_int("XDIM", 15, min = 2, max = 100),
    YDIM = env_int("YDIM", 15, min = 2, max = 100),
    TRAIN_CELLS_PER_SAMPLE = env_int("TRAIN_CELLS_PER_SAMPLE", 200000, min = 0),
    EXTRAPOLATION = env_choice("EXTRAPOLATION", "shift", c("shift", "spline")),
    MARKERS = env_list("MARKERS"),
    SEED = env_int("SEED", 42)
  )
  if (p$MODE == "controls" && !length(p$REFERENCE_VALUES)) p$REFERENCE_VALUES <- "tcu1"
  if (p$N_METACLUSTERS > p$XDIM * p$YDIM) stop("N_METACLUSTERS must not exceed XDIM * YDIM")
  p
}

norm_keyword <- function(p) {
  parts <- c("cytonorm", paste0("mode=", p$MODE))
  if (p$MODE == "controls") {
    parts <- c(parts, paste0("reference=", p$REFERENCE_COLUMN, ":", paste(p$REFERENCE_VALUES, collapse = "+")))
  }
  paste(c(parts, paste0("goal=", p$GOAL), paste0("metaclusters=", p$N_METACLUSTERS),
          paste0("quantiles=", p$N_QUANTILES), paste0("extrapolation=", p$EXTRAPOLATION),
          paste0("seed=", p$SEED)), collapse = ";")
}

seeded_rows <- function(n, k, seed) {
  if (k <= 0 || n <= k) return(seq_len(n))
  set.seed(seed)
  sort(sample.int(n, k))
}

marker_quantiles <- function(m, channels, markers, sample_id, batch, stage) {
  q <- vapply(channels, function(ch) stats::quantile(m[, ch], QC_PROBS, names = FALSE, na.rm = TRUE),
              numeric(length(QC_PROBS)))
  q <- matrix(q, nrow = length(QC_PROBS))
  data.frame(sample_id = sample_id, batch = batch,
             marker = rep(markers, each = length(QC_PROBS)),
             channel = rep(channels, each = length(QC_PROBS)),
             quantile = rep(QC_PROBS, times = length(channels)),
             stage = stage, value = as.vector(q), stringsAsFactors = FALSE)
}

#' Apply a trained CytoNorm model to one expression matrix (events x channels, colnames
#' $PnN). Returns list(exprs, metacluster, extrapolated). Same per-metacluster spline
#' application as CytoNorm::CytoNorm.normalize / QuantileNorm.normalize, infinities included.
#'
#' extrapolation = "spline": CytoNorm's own behaviour - outside the batch's outermost
#'   training quantiles the monotone spline continues linearly with its end slope. When a
#'   sample expresses a marker far beyond what the references show (e.g. CD103 on
#'   stimulated cells, while the unstimulated reference is all negative), that end slope
#'   (often 5-10) amplifies the difference and values land far off scale.
#' extrapolation = "shift" (default): outside [lowest, highest] training quantile of the
#'   batch, the value is moved by the spline's correction at that end quantile (slope 1),
#'   so the tails are shifted like their nearest trained quantile, never stretched.
#'   Identical to CytoNorm inside the trained range.
apply_cytonorm <- function(model, ff, m, label, extrapolation = "shift") {
  fsom_file <- suppressWarnings(suppressMessages(FlowSOM::NewData(model$fsom, ff)))
  mc <- FlowSOM::GetMetaclusters(fsom_file)
  n_extra <- 0
  for (cl in unique(model$fsom$metaclustering)) {
    idx <- which(as.character(mc) == as.character(cl))
    if (!length(idx)) next
    res <- model$clusterRes[[cl]]
    if (!label %in% names(res$splines)) stop(sprintf("model has no splines for batch '%s'", label))
    for (ch in res$channels) {
      spl <- res$splines[[label]][[ch]]
      v <- m[idx, ch]
      x <- spl(v)
      inf <- is.infinite(x)
      if (any(inf)) x[inf] <- sign(x[inf]) * max(abs(x[!inf]))
      if (extrapolation == "shift") {
        q <- res$quantiles[[label]][, ch]
        if (!is.null(q) && !anyNA(q)) {
          lo <- min(q)
          hi <- max(q)
          below <- v < lo
          above <- v > hi
          if (any(below)) x[below] <- v[below] + (spl(lo) - lo)
          if (any(above)) x[above] <- v[above] + (spl(hi) - hi)
          n_extra <- n_extra + sum(below) + sum(above)
        }
      }
      m[idx, ch] <- x
    }
  }
  list(exprs = m, metacluster = as.integer(as.character(mc)), extrapolated = n_extra)
}

run_cytonorm <- function(input_dir, output_dir, app_dir, p = cytonorm_params()) {
  t0 <- Sys.time()
  log_line("processor-cytonorm ", app_version(app_dir), " MODE=", p$MODE, " GOAL=", p$GOAL,
           " N_METACLUSTERS=", p$N_METACLUSTERS, " N_QUANTILES=", p$N_QUANTILES, " SEED=", p$SEED)
  suppressPackageStartupMessages({
    library(flowCore)
    library(FlowSOM)
    library(CytoNorm)
  })

  # ---- inputs ----
  fcs <- find_fcs(input_dir, "fcs")
  samples_csv <- find_named(input_dir, "samples.csv")
  panel_csv <- find_named(input_dir, "panel.csv")
  if (!length(fcs)) stop("no .fcs files found under INPUT_DIR")
  if (is.na(samples_csv)) stop("samples.csv not found under INPUT_DIR")
  if (is.na(panel_csv)) stop("panel.csv not found under INPUT_DIR")
  needed <- p$BATCH_COLUMN
  if (p$MODE == "controls") needed <- c(needed, p$REFERENCE_COLUMN)
  samples <- read_samples(samples_csv, needed)
  panel <- read_panel(panel_csv)
  sel <- select_markers(panel, p$MARKERS)
  channels <- sel$channel
  markers <- sel$marker
  clust <- sel$channel[sel$clustering]
  if (!length(clust)) clust <- channels
  log_line(sprintf("inputs: %d FCS files, samples.csv %s, panel.csv %s; %d markers normalised, %d used for clustering",
                   length(fcs), rel_path(samples_csv, input_dir), rel_path(panel_csv, input_dir),
                   length(channels), length(clust)))
  if (any(!nzchar(sel$cofactor))) {
    log_line("WARNING: panel.csv gives no cofactor for ", sum(!nzchar(sel$cofactor)),
             " selected marker(s); CytoNorm expects transformed values and applies no transform")
  }

  m <- match_samples(fcs, samples)
  m$batch <- ifelse(is.na(m$row), NA_character_, samples[[p$BATCH_COLUMN]][m$row])
  status <- list()
  for (i in which(is.na(m$row))) {
    status[[m$sample_id[i]]] <- list(status = "failed", reason = "not in samples.csv")
    log_line(m$sample_id[i], ": FAILED not in samples.csv")
  }
  missing_fcs <- setdiff(samples$sample_id, m$sample_id)
  if (length(missing_fcs)) {
    log_line(sprintf("%d sample(s) in samples.csv have no FCS file (dropped upstream?): %s",
                     length(missing_fcs), paste(missing_fcs, collapse = ", ")))
  }
  m <- m[!is.na(m$row), , drop = FALSE]
  bad_batch <- is.na(m$batch) | !nzchar(m$batch) | m$batch == "NA"
  for (i in which(bad_batch)) {
    status[[m$sample_id[i]]] <- list(status = "failed", reason = paste0("no ", p$BATCH_COLUMN))
    log_line(m$sample_id[i], ": FAILED empty ", p$BATCH_COLUMN)
  }
  m <- m[!bad_batch, , drop = FALSE]
  m <- m[order(m$batch, m$sample_id), , drop = FALSE]
  if (!nrow(m)) stop("no usable samples")

  if (p$MODE == "controls") {
    m$training <- samples[[p$REFERENCE_COLUMN]][m$row] %in% p$REFERENCE_VALUES
  } else {
    m$training <- TRUE
  }
  train_batches <- sort(unique(m$batch[m$training]))
  if (length(train_batches) < 2) {
    stop(sprintf("need training samples in at least 2 batches, found %d (MODE=%s%s)", length(train_batches), p$MODE,
                 if (p$MODE == "controls") paste0(", ", p$REFERENCE_COLUMN, " in {", paste(p$REFERENCE_VALUES, collapse = ","), "}") else ""))
  }
  no_ref <- setdiff(unique(m$batch), train_batches)
  if (length(no_ref)) log_line("WARNING: no reference sample in batch(es) ", paste(no_ref, collapse = ", "),
                               "; their samples cannot be normalised and are skipped")
  if (p$GOAL != "mean" && !p$GOAL %in% train_batches) {
    stop(sprintf("GOAL must be 'mean' or a batch with training samples (%s), got '%s'",
                 paste(train_batches, collapse = ", "), p$GOAL))
  }

  # ---- training aggregates, one FCS per batch in /tmp ----
  tmp <- file.path(tempdir(), "cytonorm")
  dir.create(file.path(tmp, "train"), recursive = TRUE, showWarnings = FALSE)
  t_train <- Sys.time()
  train_files <- character()
  train_counts <- list()
  qc_q <- list()
  qc_mc <- list()
  for (b in train_batches) {
    rows <- which(m$batch == b & m$training)
    agg <- NULL
    for (r in rows) {
      t1 <- Sys.time()
      ok <- tryCatch({
        f <- read_fcs(m$path[r])
        missing_ch <- setdiff(unique(c(channels, clust)), colnames(f$exprs))
        if (length(missing_ch)) stop("missing channel(s) ", paste(missing_ch, collapse = ", "))
        keep <- seeded_rows(nrow(f$exprs), p$TRAIN_CELLS_PER_SAMPLE, p$SEED + r)
        x <- f$exprs[keep, unique(c(channels, clust)), drop = FALSE]
        agg <- rbind(agg, x)
        train_counts[[m$sample_id[r]]] <- list(batch = b, events = nrow(f$exprs), used = length(keep))
        log_line(sprintf("train %s (%s): %d of %d events (%.1fs)", m$sample_id[r], b, length(keep),
                         nrow(f$exprs), as.numeric(difftime(Sys.time(), t1, units = "secs"))))
        TRUE
      }, error = function(e) {
        log_line(sprintf("train %s (%s): FAILED %s", m$sample_id[r], b, conditionMessage(e)))
        status[[m$sample_id[r]]] <<- list(status = "failed", reason = paste("training read failed:", conditionMessage(e)))
        FALSE
      })
    }
    if (is.null(agg) || nrow(agg) < 100) {
      log_line("WARNING: batch ", b, " has no usable training events; its samples are skipped")
      next
    }
    path <- file.path(tmp, "train", paste0(make.names(b), ".fcs"))
    suppressWarnings(flowCore::write.FCS(flowCore::flowFrame(agg), path))
    train_files[b] <- path
    rm(agg)
  }
  train_batches <- names(train_files)
  if (length(train_batches) < 2) stop("fewer than 2 batches with usable training data")
  n_train <- sum(vapply(train_counts, function(x) x$used, numeric(1)))
  log_line(sprintf("training: %d batches, %d events; FlowSOM %dx%d, %d metaclusters on %d markers",
                   length(train_batches), n_train, p$XDIM, p$YDIM, p$N_METACLUSTERS, length(clust)))

  model <- CytoNorm::CytoNorm.train(
    files = unname(train_files),
    labels = train_batches,
    channels = channels,
    transformList = NULL,
    outputDir = file.path(tmp, "model"),
    FlowSOM.params = list(nCells = min(1e6, n_train), xdim = p$XDIM, ydim = p$YDIM,
                          nClus = p$N_METACLUSTERS, scale = FALSE, colsToUse = clust),
    normMethod.train = CytoNorm::QuantileNorm.train,
    normParams = list(nQ = p$N_QUANTILES, goal = p$GOAL),
    seed = p$SEED, clean = TRUE, plot = FALSE, verbose = FALSE, recompute = TRUE,
    truncate_max_range = FALSE, transformation = FALSE)
  t_train_s <- as.numeric(difftime(Sys.time(), t_train, units = "secs"))
  log_line(sprintf("CytoNorm.train done (%.1fs)", t_train_s))

  # Metacluster share per training batch and its CV (CytoNorm's testCV check: the clustering
  # should not itself be driven by the batch effect; CV > 1.5-2 is a warning sign).
  fsom <- model$fsom
  file_col <- fsom$data[, "File"]
  train_mc <- FlowSOM::GetMetaclusters(fsom)
  tab <- table(factor(file_col, levels = seq_along(train_batches)), train_mc)
  pct <- 100 * tab / rowSums(tab)
  cv <- apply(pct, 2, function(x) stats::sd(x) / mean(x))
  cluster_cv <- data.frame(metacluster = colnames(pct), cv = round(cv, 4), stringsAsFactors = FALSE)
  for (j in seq_along(train_batches)) cluster_cv[[paste0("pct_", train_batches[j])]] <- round(as.numeric(pct[j, ]), 3)
  # Drop the aggregated training events before keeping/saving the model (NewData does not
  # need them; it maps on the SOM codes).
  model$fsom$data <- NULL
  # info$parameters holds the FlowSOM call with the training data embedded (~100 MB on the
  # demo data) and outliers is a per-training-cell table; NewData needs neither.
  model$fsom$info$parameters <- NULL
  model$fsom$outliers <- NULL
  pretty <- stats::setNames(paste0(markers, " (", channels, ")"), channels)
  model$fsom$prettyColnames[names(pretty)[names(pretty) %in% names(model$fsom$prettyColnames)]] <-
    pretty[names(pretty) %in% names(model$fsom$prettyColnames)]

  # ---- apply, one file at a time ----
  out_fcs <- file.path(output_dir, "fcs")
  out_qc <- file.path(output_dir, "qc", "cytonorm")
  dir.create(out_fcs, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_qc, recursive = TRUE, showWarnings = FALSE)
  kw_norm <- norm_keyword(p)
  counts <- list()
  t_apply <- Sys.time()
  for (r in seq_len(nrow(m))) {
    sid <- m$sample_id[r]
    b <- m$batch[r]
    if (!is.null(status[[sid]])) next
    t1 <- Sys.time()
    if (!b %in% train_batches) {
      status[[sid]] <- list(status = "failed", reason = paste0("no reference sample in batch ", b))
      log_line(sid, " (", b, "): FAILED no reference sample in its batch")
      next
    }
    res <- tryCatch({
      f <- read_fcs(m$path[r])
      missing_ch <- setdiff(unique(c(channels, clust)), colnames(f$exprs))
      if (length(missing_ch)) stop("missing channel(s) ", paste(missing_ch, collapse = ", "))
      before <- f$exprs
      a <- apply_cytonorm(model, f$flowframe, f$exprs, b, p$EXTRAPOLATION)
      f$flowframe <- NULL
      write_fcs(file.path(out_fcs, paste0(sid, ".fcs")), a$exprs, f$keywords,
                c(PENNSIEVE_NORM = kw_norm))
      changed <- colnames(before)[colSums(abs(a$exprs - before) > 0, na.rm = TRUE) > 0]
      if (length(setdiff(changed, channels))) stop("non-marker channel changed: ", paste(setdiff(changed, channels), collapse = ", "))
      counts[[sid]] <- as.data.frame(table(metacluster = a$metacluster), stringsAsFactors = FALSE)
      for (cl in sort(unique(a$metacluster))) {
        ix <- which(a$metacluster == cl)
        if (length(ix) < 50) next
        qc_mc[[length(qc_mc) + 1]] <- data.frame(
          sample_id = sid, batch = b, metacluster = cl, n = length(ix), marker = markers, channel = channels,
          before = apply(before[ix, channels, drop = FALSE], 2, stats::median),
          after = apply(a$exprs[ix, channels, drop = FALSE], 2, stats::median), stringsAsFactors = FALSE)
      }
      qc_q[[length(qc_q) + 1]] <- marker_quantiles(before, channels, markers, sid, b, "before")
      qc_q[[length(qc_q) + 1]] <- marker_quantiles(a$exprs, channels, markers, sid, b, "after")
      list(n = nrow(a$exprs), shift = mean(abs(a$exprs[, channels] - before[, channels])),
           extra = a$extrapolated / (nrow(a$exprs) * length(channels)),
           range = range(a$exprs[, channels]))
    }, error = function(e) e)
    secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    if (inherits(res, "error")) {
      status[[sid]] <- list(status = "failed", reason = conditionMessage(res))
      log_line(sprintf("%s (%s): FAILED %s (%.1fs)", sid, b, conditionMessage(res), secs))
    } else {
      status[[sid]] <- list(status = "ok", batch = b, events = res$n, reference = m$training[r],
                            mean_abs_change = round(res$shift, 5),
                            outside_trained_range = round(res$extra, 5),
                            marker_range_after = round(res$range, 3), seconds = round(secs, 2))
      log_line(sprintf("%s (%s%s): %d events normalised, mean |change| %.4f, %.2f%% of values outside the trained quantiles, marker range after [%.2f, %.2f] (%.1fs)",
                       sid, b, if (m$training[r]) ", reference" else "", res$n, res$shift, 100 * res$extra,
                       res$range[1], res$range[2], secs))
    }
    gc(verbose = FALSE)
  }
  t_apply_s <- as.numeric(difftime(Sys.time(), t_apply, units = "secs"))

  # ---- qc ----
  saveRDS(model, file.path(out_qc, "cytonorm_model.rds"))
  write_csv(cluster_cv, file.path(out_qc, "metacluster_cv.csv"))
  if (length(counts)) {
    mc_counts <- do.call(rbind, lapply(names(counts), function(s) {
      d <- counts[[s]]
      data.frame(sample_id = s, batch = m$batch[m$sample_id == s], metacluster = d$metacluster,
                 n = d$Freq, stringsAsFactors = FALSE)
    }))
    write_csv(mc_counts, file.path(out_qc, "metacluster_counts.csv"))
  }
  if (length(qc_q)) {
    q <- do.call(rbind, qc_q)
    wide <- stats::reshape(q, idvar = c("sample_id", "batch", "marker", "channel", "quantile"),
                           timevar = "stage", direction = "wide")
    names(wide) <- sub("^value\\.", "", names(wide))
    wide$reference <- wide$sample_id %in% m$sample_id[m$training]
    has_group <- p$REFERENCE_COLUMN %in% names(samples)
    wide$group <- if (has_group) samples[[p$REFERENCE_COLUMN]][match(wide$sample_id, samples$sample_id)] else NA
    wide <- wide[order(match(wide$marker, markers), wide$quantile, wide$batch, wide$sample_id), ]
    rownames(wide) <- NULL
    write_csv(wide[wide$reference, c("sample_id", "batch", "marker", "channel", "quantile", "before", "after")],
              file.path(out_qc, "reference_quantiles.csv"))
    # Spread across batches of replicate samples (same REFERENCE_COLUMN value in several
    # batches, e.g. tcu1 = training references, tcu2 = a held-out second aliquot of the same
    # donor): SD across batches of each marker quantile, before and after. Per batch,
    # replicates are averaged first. This is what CytoNorm is meant to shrink; held-out
    # replicates show whether it generalises beyond the samples it was trained on.
    groups <- if (has_group) wide else transform(wide, group = ifelse(reference, "training", NA))
    groups <- groups[!is.na(groups$group) & nzchar(groups$group), ]
    nb <- tapply(groups$batch, groups$group, function(x) length(unique(x)))
    groups <- groups[groups$group %in% names(nb)[nb >= 2], ]
    if (nrow(groups)) {
      per_batch <- stats::aggregate(cbind(before, after) ~ group + batch + marker + channel + quantile,
                                    data = groups, FUN = mean)
      spread <- stats::aggregate(cbind(before, after) ~ group + marker + channel + quantile,
                                 data = per_batch, FUN = stats::sd)
      names(spread)[match(c("before", "after"), names(spread))] <- c("sd_before", "sd_after")
      spread$n_batches <- as.integer(nb[spread$group])
      ref_groups <- unique(groups$group[groups$reference])
      spread$role <- ifelse(spread$group %in% ref_groups, "training", "held-out")
      spread$reduction <- round(1 - spread$sd_after / spread$sd_before, 4)
      spread <- spread[order(spread$role != "training", spread$group, match(spread$marker, markers), spread$quantile),
                       c("group", "role", "n_batches", "marker", "channel", "quantile", "sd_before", "sd_after", "reduction")]
      write_csv(spread, file.path(out_qc, "replicate_spread.csv"))
      if (length(qc_mc)) {
        # The same within metaclusters (what CytoNorm actually aligns): SD across batches of
        # each replicate group's per-metacluster marker median. Whole-sample quantiles of a
        # bimodal marker also move with the cell composition, which CytoNorm does not change.
        mcd <- do.call(rbind, qc_mc)
        mcd$group <- if (has_group) samples[[p$REFERENCE_COLUMN]][match(mcd$sample_id, samples$sample_id)] else
          ifelse(mcd$sample_id %in% m$sample_id[m$training], "training", NA)
        mcd <- mcd[mcd$group %in% unique(groups$group), ]
        pb <- stats::aggregate(cbind(before, after) ~ group + batch + metacluster + marker + channel, data = mcd, FUN = mean)
        nbm <- stats::aggregate(batch ~ group + metacluster + marker, data = pb, FUN = length)
        sp <- stats::aggregate(cbind(before, after) ~ group + metacluster + marker + channel, data = pb, FUN = stats::sd)
        names(sp)[match(c("before", "after"), names(sp))] <- c("sd_before", "sd_after")
        sp$n_batches <- nbm$batch[match(paste(sp$group, sp$metacluster, sp$marker), paste(nbm$group, nbm$metacluster, nbm$marker))]
        sp <- sp[sp$n_batches >= 2, ]
        sp$role <- ifelse(sp$group %in% ref_groups, "training", "held-out")
        sp$reduction <- round(1 - sp$sd_after / sp$sd_before, 4)
        sp <- sp[order(sp$role != "training", sp$group, sp$metacluster, match(sp$marker, markers)),
                 c("group", "role", "metacluster", "n_batches", "marker", "channel", "sd_before", "sd_after", "reduction")]
        write_csv(sp, file.path(out_qc, "replicate_spread_by_metacluster.csv"))
        for (g in unique(sp$group)) {
          x <- sp[sp$group == g, ]
          log_line(sprintf("replicates %s within metaclusters: SD of per-batch medians, mean over metaclusters x markers: before %.4f, after %.4f",
                           g, mean(x$sd_before), mean(x$sd_after)))
        }
      }
      for (g in unique(spread$group)) {
        med <- spread[spread$group == g & spread$quantile == 0.5, ]
        log_line(sprintf("replicates %s (%s, %d batches): SD of per-batch medians across batches, mean over markers: before %.4f, after %.4f",
                         g, med$role[1], med$n_batches[1], mean(med$sd_before), mean(med$sd_after)))
      }
    }
  }
  dir.create(file.path(out_qc, "splines"), showWarnings = FALSE)
  for (b in train_batches) {
    ok <- tryCatch({
      pl <- CytoNorm::plotSplines(model, batches = b, channels = channels, groupClusters = TRUE)
      ncol <- ceiling(sqrt(length(channels)))
      ggplot2::ggsave(file.path(out_qc, "splines", paste0(make.names(b), ".png")), pl[[1]],
                      width = 2.2 * ncol, height = 2.2 * ceiling(length(channels) / ncol), dpi = 90,
                      limitsize = FALSE)
      TRUE
    }, error = function(e) { log_line("WARNING: splines plot for ", b, " failed: ", conditionMessage(e)); FALSE })
  }

  ok_n <- sum(vapply(status, function(s) identical(s$status, "ok"), logical(1)))
  failed_n <- length(status) - ok_n
  prov <- list(
    app = "processor-cytonorm", version = app_version(app_dir),
    tool = "CytoNorm", tool_versions = package_versions(CYTONORM_PKGS),
    cytonorm_source = "https://github.com/saeyslab/CytoNorm@278f352d8edffd7a5d3e6b0b76c402d991a1e7a9",
    parameters = p,
    keyword = list(PENNSIEVE_NORM = kw_norm),
    normalised_markers = data.frame(channel = channels, marker = markers, clustering = channels %in% clust),
    training = list(batches = train_batches, samples = train_counts, flowsom_cells = min(1e6, n_train)),
    inputs = list(fcs = rel_path(fcs, input_dir), samples = rel_path(samples_csv, input_dir),
                  panel = rel_path(panel_csv, input_dir)),
    timings = list(started = format(t0, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), finished = iso_now(),
                   train_seconds = round(t_train_s, 1), apply_seconds = round(t_apply_s, 1),
                   total_seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)),
    samples = status, n_ok = ok_n, n_failed = failed_n
  )
  write_json(prov, file.path(out_fcs, "provenance.json"))
  write_json(prov, file.path(out_qc, "provenance.json"))
  log_line(sprintf("done: %d sample(s) normalised, %d failed, %.1fs", ok_n, failed_n,
                   as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (ok_n == 0) stop("every sample failed")
  invisible(prov)
}

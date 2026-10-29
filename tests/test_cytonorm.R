# Unit tests on a tiny synthetic fixture (3 batches x 3 samples x 3000 events).

TEST_ENV <- c(XDIM = "5", YDIM = "5", N_METACLUSTERS = "3", N_QUANTILES = "51",
              REFERENCE_VALUES = "ctrl", SEED = "7")

test_fcs_roundtrip <- function() {
  d <- new_tmp("rt")
  m <- cbind(`FSC-A` = c(1, 2.5, 1e6), `PE-A` = c(-1.25, 0, 3.5))
  kw <- c("$P1N" = "FSC-A", "$P2N" = "PE-A", "$P2S" = "CD19", "$P1R" = "262144", "$P2R" = "262144",
          "CUSTOM|KEY" = "a|b", "$DATE" = "01-Jan-2024", "$BEGINDATA" = "999")
  write_fcs(file.path(d, "x.fcs"), m, kw, c(PENNSIEVE_NORM = "test"))
  hdr <- readChar(file.path(d, "x.fcs"), 6, useBytes = TRUE)
  expect_true(hdr == "FCS3.1", "FCS 3.1 header")
  f <- read_fcs(file.path(d, "x.fcs"))
  expect_equal_num(f$exprs, m, 1e-6 * 1e6, "values round-trip (float32)")
  expect_true(identical(unname(f$markers), c("", "CD19")), "PnS kept")
  expect_true(fcs_keyword(f$keywords, "CUSTOM|KEY") == "a|b", "delimiter in keyword escaped")
  expect_true(fcs_keyword(f$keywords, "PENNSIEVE_NORM") == "test", "extra keyword written")
  expect_true(fcs_keyword(f$keywords, "$DATATYPE") == "F", "float32")
  expect_true(fcs_keyword(f$keywords, "$P2B") == "32", "PnB 32")
  expect_true(fcs_keyword(f$keywords, "$TOT") == "3", "TOT")
  expect_error(write_fcs(file.path(d, "y.fcs"), m[, 2:1], kw), "do not match")
}

test_controls_mode <- function() {
  root <- new_tmp("ctl")
  inp <- file.path(root, "in")
  out <- file.path(root, "out")
  make_synthetic(inp)
  prov <- with_env(TEST_ENV, run_cytonorm(inp, out, "R", cytonorm_params()))
  expect_true(prov$n_ok == 9 && prov$n_failed == 0, "all 9 samples ok")
  for (s in c("b1_ctrl", "b2_s1", "b3_s2")) {
    a <- read_fcs(file.path(inp, "fcs", paste0(s, ".fcs")))
    b <- read_fcs(file.path(out, "fcs", paste0(s, ".fcs")))
    expect_true(identical(colnames(a$exprs), colnames(b$exprs)), "same channels")
    expect_true(identical(a$markers, b$markers), "same PnS")
    expect_true(nrow(a$exprs) == nrow(b$exprs), "no events dropped")
    for (ch in c("FSC-A", "SSC-A", "AF-A", "Time")) expect_equal_num(a$exprs[, ch], b$exprs[, ch], 0, ch)
    expect_true(fcs_keyword(b$keywords, "GROUPNAME") == "donor A", "keywords kept")
    expect_true(fcs_keyword(b$keywords, "PENNSIEVE_PREP") == "transform=arcsinh;cofactor=6000", "prep keyword kept")
    expect_true(grepl("^cytonorm;mode=controls;reference=control:ctrl;goal=mean", fcs_keyword(b$keywords, "PENNSIEVE_NORM")),
                "PENNSIEVE_NORM keyword")
  }
  spread <- read.csv(file.path(out, "qc", "cytonorm", "replicate_spread.csv"))
  ctrl <- spread[spread$group == "ctrl" & spread$quantile == 0.5, ]
  expect_true(mean(ctrl$sd_after) < 0.25 * mean(ctrl$sd_before), "reference medians aligned across batches")
  held <- spread[spread$group == "s1" & spread$quantile == 0.5, ]
  expect_true(all(held$role == "held-out") && mean(held$sd_after) < 0.5 * mean(held$sd_before),
              "held-out replicate s1 aligned too")
  for (f in c("cytonorm_model.rds", "replicate_spread_by_metacluster.csv", "metacluster_cv.csv", "metacluster_counts.csv", "reference_quantiles.csv",
              "provenance.json", "splines/b1.png")) {
    expect_true(file.exists(file.path(out, "qc", "cytonorm", f)), paste("qc file", f))
  }
  expect_true(file.exists(file.path(out, "fcs", "provenance.json")), "fcs provenance")
  rq <- read.csv(file.path(out, "qc", "cytonorm", "reference_quantiles.csv"))
  expect_true(all(grepl("_ctrl$", rq$sample_id)), "reference_quantiles only has references")
}

test_matches_cytonorm_normalize <- function() {
  # The in-memory application must give what CytoNorm.normalize writes.
  root <- new_tmp("eq")
  inp <- file.path(root, "in")
  out <- file.path(root, "out")
  make_synthetic(inp)
  with_env(c(TEST_ENV, EXTRAPOLATION = "spline"), run_cytonorm(inp, out, "R", cytonorm_params()))
  model <- readRDS(file.path(out, "qc", "cytonorm", "cytonorm_model.rds"))
  f <- file.path(inp, "fcs", "b2_s1.fcs")
  ref_dir <- file.path(root, "ref")
  suppressWarnings(suppressMessages(CytoNorm::CytoNorm.normalize(
    model = model, files = f, labels = "b2", transformList = NULL, transformList.reverse = NULL,
    outputDir = ref_dir, prefix = "Norm_", clean = TRUE, verbose = FALSE,
    truncate_max_range = FALSE, transformation = FALSE)))
  ref <- flowCore::exprs(flowCore::read.FCS(file.path(ref_dir, "Norm_b2_s1.fcs"), transformation = FALSE,
                                            truncate_max_range = FALSE))
  ours <- read_fcs(file.path(out, "fcs", "b2_s1.fcs"))$exprs
  for (ch in c("BV421-A", "PE-A", "APC-A", "FITC-A")) expect_equal_num(ours[, ch], ref[, ch], 1e-4, ch)
}

test_cohort_mode_and_goal <- function() {
  root <- new_tmp("coh")
  inp <- file.path(root, "in")
  make_synthetic(inp)
  env <- TEST_ENV
  env["MODE"] <- "cohort"
  env["GOAL"] <- "b1"
  prov <- with_env(env, run_cytonorm(inp, file.path(root, "out"), "R", cytonorm_params()))
  expect_true(prov$n_ok == 9, "cohort: all samples ok")
  expect_true(length(prov$training$samples) == 9, "cohort: every sample trains")
  b1_in <- read_fcs(file.path(inp, "fcs", "b1_s1.fcs"))$exprs
  b1_out <- read_fcs(file.path(root, "out", "fcs", "b1_s1.fcs"))$exprs
  expect_true(mean(abs(b1_in[, "BV421-A"] - b1_out[, "BV421-A"])) < 0.05, "GOAL=b1 leaves b1 almost unchanged")
  env["GOAL"] <- "nope"
  expect_error(with_env(env, run_cytonorm(inp, file.path(root, "out2"), "R", cytonorm_params())), "GOAL")
}

test_partial_failures_and_unlinked_layout <- function() {
  root <- new_tmp("fail")
  src <- file.path(root, "src")
  make_synthetic(src, extra_sample = "stray", missing_ref_batch = "b3")
  # Unlinked layout: everything side by side at the INPUT_DIR root, as symlinks.
  inp <- file.path(root, "in")
  dir.create(inp)
  for (f in list.files(file.path(src, "fcs"), full.names = TRUE)) file.symlink(f, file.path(inp, basename(f)))
  for (f in list.files(file.path(src, "metadata"), full.names = TRUE)) file.symlink(f, file.path(inp, basename(f)))
  prov <- with_env(TEST_ENV, run_cytonorm(inp, file.path(root, "out"), "R", cytonorm_params()))
  expect_true(prov$samples$stray$status == "failed", "sample missing from samples.csv fails")
  expect_true(prov$samples$b3_s1$status == "failed", "batch without reference fails")
  expect_true(prov$samples$b1_s1$status == "ok", "others ok")
  expect_true(!file.exists(file.path(root, "out", "fcs", "b3_s1.fcs")), "no output for failed sample")
  expect_error(with_env(c(TEST_ENV, REFERENCE_VALUES = "none"),
                        run_cytonorm(inp, file.path(root, "out3"), "R", cytonorm_params())), "at least 2 batches")
}

test_params <- function() {
  expect_error(with_env(c(MODE = "bogus"), cytonorm_params()), "MODE")
  expect_error(with_env(c(N_QUANTILES = "x"), cytonorm_params()), "integer")
  p <- with_env(c(REFERENCE_VALUES = "a, b", MARKERS = "CD3,CD19"), cytonorm_params())
  expect_true(identical(p$REFERENCE_VALUES, c("a", "b")) && identical(p$MARKERS, c("CD3", "CD19")), "lists parsed")
  panel <- data.frame(channel = c("A", "B"), marker = c("CD3", "B"), type = c("marker", "scatter"),
                      cofactor = "", clustering = TRUE, state = FALSE)
  expect_true(identical(select_markers(panel)$channel, "A"), "default markers = type marker")
  expect_error(select_markers(panel, "CD99"), "not in panel")
}

test_extrapolation_shift <- function() {
  # A sample with marker values far above anything the references show: CytoNorm's spline
  # extrapolates with its end slope; EXTRAPOLATION=shift moves those values by the end
  # quantile's correction instead.
  root <- new_tmp("ext")
  inp <- file.path(root, "in")
  make_synthetic(inp)
  f <- read_fcs(file.path(inp, "fcs", "b2_s2.fcs"))
  x <- f$exprs
  hot <- seq(1, nrow(x), by = 10)
  x[hot, "FITC-A"] <- x[hot, "FITC-A"] + 6  # GZMB far beyond the references' 99th percentile
  write_fcs(file.path(inp, "fcs", "b2_s2.fcs"), x, f$keywords)
  for (mode in c("shift", "spline")) {
    with_env(c(TEST_ENV, EXTRAPOLATION = mode), run_cytonorm(inp, file.path(root, mode), "R", cytonorm_params()))
  }
  sh <- read_fcs(file.path(root, "shift", "fcs", "b2_s2.fcs"))$exprs[, "FITC-A"]
  sp <- read_fcs(file.path(root, "spline", "fcs", "b2_s2.fcs"))$exprs[, "FITC-A"]
  ch_shift <- sh[hot] - x[hot, "FITC-A"]
  expect_true(max(abs(ch_shift)) < 1, "shift: out-of-range values moved by an end-quantile offset only")
  inside <- read_fcs(file.path(root, "shift", "fcs", "b1_ctrl.fcs"))$exprs[, "PE-A"]
  inside_sp <- read_fcs(file.path(root, "spline", "fcs", "b1_ctrl.fcs"))$exprs[, "PE-A"]
  # Only events beyond their metacluster's outermost training quantiles (~1 % per tail) differ.
  expect_true(mean(abs(inside - inside_sp) > 1e-5) < 0.05, "shift = CytoNorm inside the trained range")
  expect_true(grepl("extrapolation=shift", fcs_keyword(read_fcs(file.path(root, "shift", "fcs", "b1_ctrl.fcs"))$keywords, "PENNSIEVE_NORM")),
              "keyword records extrapolation")
}

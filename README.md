# processor-cytonorm

Batch normalisation of flow / mass cytometry data with **CytoNorm** (Van Gassen et al. 2020;
CytoNorm 2.0, Quintelier et al. 2025), one of the two correction branches of the
*Multi-batch Harmonization & QC* workflow:

```
fcs-prep -> peacoqc -> { batch-assess(raw), cytonorm -> batch-assess(cytonorm), cycombine -> batch-assess(cycombine) } -> cytometry-report
```

**Why.** When a study is stained and acquired over several days (batches), each batch
shifts and stretches the marker intensities a little differently: instrument settings, a
new antibody lot, staining time. Those technical differences can be as large as the
biology and show up as fake differences between patients measured on different days.
CytoNorm removes them using a sample that should look the same every time - the same
control aliquot (here donor 181C, unstimulated) measured in every batch. Whatever differs
between the copies of the control is technical, and is undone in every sample of that
batch.

**How.** The batch effect is often cell-type specific (a dim marker on T cells may move
differently from a bright one on B cells), so CytoNorm first groups cells into a few
broad FlowSOM metaclusters. Then, per metacluster, batch and marker it computes 101
quantiles of the control's distribution and fits a monotone spline that maps the batch's
quantiles onto the goal quantiles (the mean over batches, or one chosen batch). Every
cell of every file is assigned to its metacluster and passed through the spline of its
batch. Cluster shares are not changed; only marker values move. With `MODE=cohort`
(CytoNorm 2.0) there is no dedicated control: each batch's pooled samples serve as its
"control", which is valid when every batch holds a comparable mix of samples.

The input values are already arcsinh-transformed by `fcs-prep`, so CytoNorm's transform
lists are identity; scatter, `Time` and other non-marker channels are written back
bit-for-bit, as are all keywords.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `fcs` | input (subfolder) | `<sample_id>.fcs` | Arcsinh-transformed FCS from `fcs-prep` / `peacoqc`. Looked for in `INPUT_DIR/fcs/`, else anywhere under `INPUT_DIR`. |
| `metadata` | input (subfolder) | `samples.csv`, `panel.csv` | `samples.csv`: `sample_id`, `BATCH_COLUMN`, `REFERENCE_COLUMN` (MODE=controls). `panel.csv`: `channel`, `marker`, `type`, `clustering`. Looked for in `INPUT_DIR/metadata/`, else by name anywhere. |
| `fcs` | output | `fcs/<sample_id>.fcs`, `fcs/provenance.json` | Same names, events, channels and keywords as the input (FCS 3.1, float32); marker channels normalised; keyword `PENNSIEVE_NORM` added, e.g. `cytonorm;mode=controls;reference=control:tcu1;goal=mean;metaclusters=5;quantiles=101;extrapolation=shift;seed=42`. |
| `qc` | output | `qc/cytonorm/...` | Method diagnostics, the model and `provenance.json`, all under a `cytonorm/` folder: in workflow B the report links the qc outputs of cytonorm and cycombine into one input port, where entries with the same name (`provenance.json`) would collide. |

### `qc/cytonorm/` files

| File | Columns / content |
|---|---|
| `reference_quantiles.csv` | `sample_id, batch, marker, channel, quantile, before, after`: the 1/5/25/50/75/95/99 % quantiles of each marker in each reference (training) sample, whole sample, before and after normalisation. |
| `replicate_spread.csv` | `group, role, n_batches, marker, channel, quantile, sd_before, sd_after, reduction`: for every value of `REFERENCE_COLUMN` that occurs in 2+ batches, the SD across batches of the per-batch quantile (replicates in a batch averaged). `role` = `training` (the references) or `held-out` (e.g. `tcu2`, a second aliquot of the same donor that the model never saw - the honest check). `reduction` = 1 - sd_after / sd_before. |
| `replicate_spread_by_metacluster.csv` | `group, role, metacluster, n_batches, marker, channel, sd_before, sd_after, reduction`: the same within each metacluster, on per-sample marker medians - what CytoNorm actually aligns. Whole-sample quantiles of bimodal markers (CD3, CD4) also move with the cell composition, which normalisation does not (and should not) change. |
| `metacluster_cv.csv` | `metacluster, cv, pct_<batch>...`: share of each metacluster in each batch's training cells and its coefficient of variation. CytoNorm assumes the clustering is not driven by the batch effect; the authors treat CVs above ~1.5-2 as a warning. |
| `metacluster_counts.csv` | `sample_id, batch, metacluster, n`: cells per metacluster in every output file. |
| `splines/<batch>.png` | CytoNorm's `plotSplines`: per marker, the batch's quantiles (x) against the goal quantiles (y) and the fitted spline, one colour per metacluster. A straight diagonal = no correction. |
| `cytonorm_model.rds` | The trained model (`CytoNorm.train` result; FlowSOM training events, outlier table and stored call removed; 14 MB on the demo data). `readRDS()` it and pass it to `CytoNorm::CytoNorm.normalize` to apply the same normalisation to new files of the same batches. |
| `provenance.json` | App and package versions, parameters, inputs, training samples and event counts, timings, per-sample status and mean absolute change. |

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `MODE` | `controls` \| `cohort` | `controls` | `controls`: train on the reference samples (CytoNorm 1; a reference is needed in every batch, batches without one are skipped and their samples fail). `cohort`: train on all samples of each batch pooled (CytoNorm 2.0, no dedicated controls). |
| `BATCH_COLUMN` | string | `batch` | `samples.csv` column with each sample's batch. |
| `REFERENCE_COLUMN` | string | `control` | `samples.csv` column that marks the references. Other values of it that recur across batches are reported as held-out replicates. |
| `REFERENCE_VALUES` | string | `tcu1` | Comma-separated values of `REFERENCE_COLUMN` that are references. Several references in one batch are pooled. |
| `GOAL` | string | `mean` | Goal distribution: `mean` (per quantile, mean over batches) or a batch id (that batch is left as is; the others are mapped onto it). CytoNorm also accepts fixed quantile values; not exposed. |
| `N_METACLUSTERS` | integer | `5` | FlowSOM metaclusters with their own splines (3-100). Few and broad: rare clusters give too few cells for 101 quantiles (CytoNorm needs 50 per cluster and batch, else the identity is used). |
| `N_QUANTILES` | integer | `101` | Quantiles per metacluster and marker (`nQ`); the splines go through them. |
| `XDIM`, `YDIM` | integer | `15`, `15` | FlowSOM grid before metaclustering (CytoNorm's default). |
| `TRAIN_CELLS_PER_SAMPLE` | integer | `200000` | At most this many events per training sample (seeded) build the model; 0 = all. Every event is still normalised. |
| `EXTRAPOLATION` | `shift` \| `spline` | `shift` | Values beyond the outermost trained quantile (1 % / 99 %) of their batch and metacluster. `spline` = CytoNorm's own behaviour: the spline continues linearly with its end slope. `shift`: such values are moved by the spline's correction at that end quantile (slope 1). See below. |
| `MARKERS` | string | (all `type=marker`) | Comma-separated markers or channels to normalise. Clustering uses those with `clustering=true` in `panel.csv` (all, if none). |
| `SEED` | integer | `42` | Seed for the training subsample and FlowSOM. |

## The demo data and its reference sample

On the CytoBatchFlagR dataset (6 batches) the references are `control = tcu1`. The Zenodo
description says the unstimulated technical controls tcu1 and tcu2 were aliquoted from
the same blood-donor leukocytes and frozen, then thawed and stained with every batch. The
FCS keywords confirm it is one donor throughout: `$TUBENAME` is `181C` for every tcu1 file
and `181C 2` / `181C  2` for every tcu2 file, in all six batches (`$GROUPNAME` is `181C` in
batch 01, `Blood Controls` elsewhere). So tcu1 is a replicate of the same sample in every
batch - what CytoNorm 1 needs - and tcu2, a second aliquot of the same donor, is a
held-out validation sample the model never sees. The stimulated controls are **not**
replicates: tcs3/tcs4 are donor `199C` in batches 01 and 05 and donor `198C` in batches
02, 03, 04 and 06 (`$TUBENAME` "199C d7 stim" / "198C d7 stim"), so their differences
between batches are partly biological.

## Values outside the trained range (`EXTRAPOLATION`)

CytoNorm's splines (`stats::splinefun(method = "monoH.FC")`) are fitted through the 1-99 %
quantiles of each batch's references and continue linearly with the end slope beyond
them. When a sample expresses a marker far beyond anything the reference shows, that end
slope multiplies the excess. On the demo data the unstimulated reference is CD103
negative (99th percentile 0.19 on the arcsinh scale in batch 01) while the stimulated
samples reach 4.5; the end slope of ~10 sent 4,183 CD103 values of batch01_tcs3 to 20-43
(arcsinh scale), clearly artefacts. CytoNorm offers `limit` (extra fixed points mapped to
themselves) but it does not bound the extrapolation beyond those points. The app
therefore applies, by default, a slope-1 continuation: a value beyond the trained range is
moved by the same amount as the end quantile itself (`EXTRAPOLATION=shift`). Inside the
trained range the result is identical to CytoNorm (a unit test checks this against
`CytoNorm.normalize`); with `EXTRAPOLATION=spline` it is identical everywhere. The share
of values outside the trained range is logged and stored per sample in
`provenance.json` (`outside_trained_range`): ~2 % for the references themselves (the
1 %/99 % tails), 6-22 % for the stimulated samples on the demo data - a sign that
references unlike the samples limit what CytoNorm can learn for those markers.

## Results on the demo data

Full demo data from `fcs-prep` (24 files, 10,212,112 events; 25 markers normalised, 14 of
them used for clustering per `panel.csv`; `data/work/prep`, no PeacoQC yet), default parameters (`MODE=controls`,
`REFERENCE_VALUES=tcu1`, `GOAL=mean`, 5 metaclusters, 101 quantiles, 200,000 training
events per reference). All 24 files normalised.

Spread of the replicate samples across the 6 batches (SD across batches of the per-batch
marker median within each metacluster, mean over 5 metaclusters x 25 markers;
`replicate_spread_by_metacluster.csv`):

| Replicate | Role | Before | After | Reduction |
|---|---|---|---|---|
| tcu1 (donor 181C) | training reference | 0.208 | 0.033 | 84 % |
| tcu2 (donor 181C, 2nd aliquot) | held out | 0.208 | 0.106 | 49 % |
| tcs3 (stimulated; donor 199C or 198C) | not a true replicate | 0.183 | 0.200 | -9 % |
| tcs4 (stimulated; donor 199C or 198C) | not a true replicate | 0.183 | 0.199 | -9 % |

The held-out tcu2 aliquot is the honest measure: about half of its between-batch spread is
removed. Per marker (tcu2, before -> after): CD45RA 0.34 -> 0.22, CD8 0.38 -> 0.14,
CD14 0.35 -> 0.11, GZMB 0.42 -> 0.14, HLADR 0.26 -> 0.12, CD16 0.31 -> 0.12; CD3, CCR6
and CD45 do not improve. The stimulated samples get no closer, as expected: they are
different donors in different batches (see above) and were not trained on; their
markers often lie outside the range of the unstimulated references.

Whole-sample medians (`replicate_spread.csv`, quantile 0.5): tcu1 0.179 -> 0.077, tcu2
0.173 -> 0.120. Metacluster shares are not batch-driven (`metacluster_cv.csv`: CV
0.27-0.48 for the four populated metaclusters; metacluster 5 holds < 0.06 % of cells and
gets the identity in most batches). Batch 01's reference has a visibly different
composition (metacluster 4: 36 % vs 16-23 % elsewhere), plausibly the separate
acquisition group (`$GROUPNAME` 181C) - CytoNorm leaves composition alone. Output values
stay in [-6.1, 7.8] on the arcsinh scale. Full outputs:
`workflow-immune-cytometry/data/work/cytonorm/`.

## Runtime

`computeTypes: [standard]`, 2 vCPU / 8 GiB, 2 h timeout. R is single-threaded here.

Measured on the full demo data above (24 files, 10.21 M events, 35 channels, 1.3 GB),
linux/amd64 under emulation on an Apple-silicon Mac (Docker Desktop, `--cpus 4`), cgroup
`memory.peak` (includes page cache of the files read):

| Stage | Wall time |
|---|---|
| Read 6 references, 1.06 M training events | 5 s |
| `CytoNorm.train` (FlowSOM 15x15 on 1.06 M cells, splits, splines) | 123 s |
| Normalise + write 24 files (2-24 s per file, ~40 s per million events) | 430 s |
| QC tables and plots | 15 s |
| **Total** | **580 s**, peak memory **2.7-3.8 GB** (two runs) |

Native amd64 should be several times faster. Memory is dominated by the largest single
file (600 k events x 35 channels here: the file, its mapping to FlowSOM and a copy) plus
training (`TRAIN_CELLS_PER_SAMPLE` x references, FlowSOM capped at 1 M cells); it does not
grow with the number of files. 8 GiB leaves room for files up to roughly 2-3 M events.
The 20,000-events-per-file fixture runs in ~80 s with a 0.7 GB peak.

The image is 600 MB (R 4.6.1, flowCore, FlowSOM, CytoNorm and their dependencies).

## Limitations

- CytoNorm corrects intensities, not cell counts: a batch effect that changes which cells
  pass gating or QC is not corrected.
- It assumes the FlowSOM metaclusters are the same cell types in every batch. Check
  `metacluster_cv.csv`; if a metacluster's share swings between batches, lower
  `N_METACLUSTERS` or restrict clustering to robust lineage markers in `panel.csv`.
- `MODE=controls` cannot normalise a batch that has no reference sample (its samples fail
  and are listed in `provenance.json`). `MODE=cohort` needs balanced batches; on the demo
  data the stimulated samples come from two different donors depending on the batch, so
  cohort mode would partly treat donor differences as batch effect.
- The splines are fitted between the 1 % and 99 % quantiles (CytoNorm 2: the 0 and 100 %
  quantiles are no longer used); values beyond them are extrapolated linearly by the
  spline, which can stretch outliers. `TRAIN_CELLS_PER_SAMPLE` subsamples large references.
- Only data on an already transformed scale make sense here (`fcs-prep` output); raw
  linear intensities are not transformed by this app.

## Running it locally

```sh
pip install -r requirements-dev.txt
make fixture DEMO=../workflow-immune-cytometry/data/zenodo-15388817   # 20,000 events per file
make test                                                             # unit tests (synthetic data)
make run                                                              # docker compose over data/fixture
```

Or by hand:

```sh
docker build --platform linux/amd64 -t processor-cytonorm .
docker run --rm --platform linux/amd64 \
  -e MODE=controls -e REFERENCE_VALUES=tcu1 \
  -v "$PWD/data/fixture:/data/input:ro" -v "$PWD/data/output:/data/output" \
  processor-cytonorm
```

`make test` runs `tests/run_tests.R` with a local R that has the packages (flowCore,
FlowSOM, CytoNorm, ggplot2, jsonlite); without one it runs the same tests in the image.
The tests build a synthetic 3-batch fixture and check that the output keeps channels,
keywords and non-marker values exactly, that replicate medians align, that the in-memory
application equals `CytoNorm.normalize`, cohort mode, `GOAL`, unlinked inputs and per-sample
failures.

## Licences

This app's code is Apache-2.0 (`LICENSE`). The image bundles R and R packages under their
own licences, several of them GPL (R, CytoNorm, FlowSOM) and flowCore's backend cytolib
AGPL-3.0, so the image as a whole is distributed under GPL/AGPL terms; see `THIRD_PARTY.md`.

## Citations

- Van Gassen S, Gaudilliere B, Angst MS, Saeys Y, Aghaeepour N. CytoNorm: A Normalization
  Algorithm for Cytometry Data. *Cytometry A* 97(3):268-278 (2020).
  doi:10.1002/cyto.a.23904. PMID 31633883.
- Quintelier KLA, Willemsen M, Bosteels V, Aerts JGJV, Saeys Y, Van Gassen S. CytoNorm 2.0:
  A flexible normalization framework for cytometry data without requiring dedicated
  controls. *Cytometry A* 107(2):69-87 (2025). doi:10.1002/cyto.a.24910. PMID 39871681.
- Van Gassen S, Callebaut B, Van Helden MJ, Lambrecht BN, Demeester P, Dhaene T, Saeys Y.
  FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry
  data. *Cytometry A* 87(7):636-645 (2015). doi:10.1002/cyto.a.22625. PMID 25573116.
- Hahne F, LeMeur N, Brinkman RR, Ellis B, Haaland P, Sarkar D, Spidlen J, Strain E,
  Gentleman R. flowCore: a Bioconductor package for high throughput flow cytometry.
  *BMC Bioinformatics* 10:106 (2009). doi:10.1186/1471-2105-10-106. PMID 19358741.
- Demo data: Koenig Z, Andorf S, Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow
  cytometry dataset. Zenodo (2025). doi:10.5281/zenodo.15388817 (CC BY 4.0).

Software: CytoNorm https://github.com/saeyslab/CytoNorm at commit `278f352d` (version 2.0.12).

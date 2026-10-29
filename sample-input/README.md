# Sample input

No FCS files are committed. The app expects the output of the immune cytometry
workflow's `fcs-prep` step (optionally after `peacoqc`): two input folders,

```
INPUT_DIR/fcs/<sample_id>.fcs           arcsinh-transformed FCS 3.1 files, one per sample
INPUT_DIR/metadata/samples.csv          sample_id, file, batch, condition, n_events (+ control, ...)
INPUT_DIR/metadata/panel.csv            channel, marker, type, cofactor, clustering, state
```

(see `workflow-immune-cytometry/docs/conventions.md`). Without links, everything may also
sit side by side at the root of INPUT_DIR; files are then found by name.

## Demo data

The CytoBatchFlagR spectral flow dataset (Koenig Z, Andorf S, Tilburgs T, Eswar S. Zenodo 2025,
doi:10.5281/zenodo.15388817, CC BY 4.0): 6 batches x 4 control samples (tcu1, tcu2 =
unstimulated aliquots of one donor; tcs3, tcs4 = anti-CD3/CD28 + IL-2 stimulated), 24 files,
10.4 M events. Download the 24 FCS files from https://zenodo.org/records/15388817 into a
folder together with a `samples.csv` sample sheet (`file, sample_id, batch, condition,
control`; one row per file, e.g. `B01_tcu1_CD45+CD14-.fcs,batch01_tcu1,batch01,unstimulated,tcu1`).

Then either run `fcs-prep` on it, or build a stand-in with the fixture script, which
subsamples each file, drops channels missing from some files (BV785-A / CD28 is absent in
batches 05-06) and arcsinh-transforms the marker channels with cofactor 6000:

```sh
pip install -r requirements-dev.txt
make fixture DEMO=/path/to/zenodo-15388817 EVENTS=20000   # -> data/fixture/{fcs,metadata}
make run                                                  # -> data/output/{fcs,qc}
```

`EVENTS=0` keeps every event (1.3 GB).

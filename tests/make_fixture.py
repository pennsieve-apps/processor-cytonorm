#!/usr/bin/env python3
"""
Build an input fixture that follows the immune cytometry file contract
(workflow-immune-cytometry/docs/conventions.md) from the raw CytoBatchFlagR demo files,
standing in for the output of fcs-prep -> peacoqc until those steps exist:

  <out>/fcs/<sample_id>.fcs        FCS 3.1, float32, marker channels arcsinh(x / 6000),
                                   scatter/Time/other channels as acquired, channels
                                   missing from any sample dropped (CD28 / BV785-A)
  <out>/metadata/samples.csv       sample_id, file, batch, condition, n_events, + extra
                                   columns of the sample sheet (control)
  <out>/metadata/panel.csv         channel, marker, type, cofactor, clustering, state

Events are subsampled per file (seeded, original order kept); --events 0 keeps all.
QC (PeacoQC) is not applied: the fixture is "prep" output, which is what peacoqc would
pass through for a clean file.

Usage: python3 tests/make_fixture.py --source <zenodo dir> --out data/fixture [--events 20000]
Needs numpy and flowio (requirements-dev.txt).
"""
from __future__ import annotations

import argparse
import csv
import os
import re
from array import array
from pathlib import Path

import flowio
import numpy as np

COFACTOR = 6000.0
# Lineage markers (used for clustering by default); everything else that is a marker is
# a functional / state marker. Matches what fcs-prep's PANEL_OVERRIDES defaults are meant to do.
LINEAGE = {
    "CD3", "CD4", "CD8", "CD19", "CD56", "CD16", "CD14", "CD45", "CD45RA", "CCR7", "CD27",
    "HLADR", "CD95", "CCR6", "CD103", "CD163",
}
SEGMENT_KEYS = {
    "$BEGINANALYSIS", "$ENDANALYSIS", "$BEGINSTEXT", "$ENDSTEXT", "$BEGINDATA", "$ENDDATA",
    "$NEXTDATA", "$PAR", "$TOT", "$DATATYPE", "$BYTEORD", "$MODE",
}


def channel_type(name: str, marker: str) -> str:
    if name.lower() == "time":
        return "time"
    if re.match(r"^(FSC|SSC)", name, re.I):
        return "scatter"
    if marker and marker != name:
        return "marker"
    return "other"


def read_fcs(path: Path):
    fd = flowio.FlowData(str(path))
    n_par = int(fd.text["par"])
    names = [fd.text[f"p{i}n"] for i in range(1, n_par + 1)]
    markers = [fd.text.get(f"p{i}s", "").strip() for i in range(1, n_par + 1)]
    data = np.asarray(fd.as_array(preprocess=False), dtype=np.float64)
    return fd, names, markers, data


def sanitise(stem: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]+", "_", stem).strip("_")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True, type=Path, help="folder with the demo FCS files + samples.csv")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--events", type=int, default=20000, help="events per file (0 = all)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--only", default="", help="comma-separated sample_ids to include (default all)")
    args = ap.parse_args()

    with open(args.source / "samples.csv", newline="") as fh:
        sheet = list(csv.DictReader(fh))
    only = {s for s in args.only.split(",") if s}
    if only:
        sheet = [r for r in sheet if r.get("sample_id") in only]

    # Pass 1: channel intersection (by $PnN) in the first file's order.
    headers = {}
    for row in sheet:
        fd = flowio.FlowData(str(args.source / row["file"]), only_text=True)
        n_par = int(fd.text["par"])
        headers[row["file"]] = [fd.text[f"p{i}n"] for i in range(1, n_par + 1)]
    first = headers[sheet[0]["file"]]
    common = [c for c in first if all(c in h for h in headers.values())]

    (args.out / "fcs").mkdir(parents=True, exist_ok=True)
    (args.out / "metadata").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    extra_cols = [c for c in sheet[0].keys() if c not in ("file", "sample_id", "batch", "condition")]
    out_rows = []
    panel = None
    for row in sheet:
        src = args.source / row["file"]
        fd, names, markers, data = read_fcs(src)
        idx = [names.index(c) for c in common]
        data = data[:, idx]
        mk = [markers[i] for i in idx]
        if args.events and data.shape[0] > args.events:
            keep = np.sort(rng.choice(data.shape[0], args.events, replace=False))
            data = data[keep]
        types = [channel_type(n, m) for n, m in zip(common, mk)]
        for j, t in enumerate(types):
            if t == "marker":
                data[:, j] = np.arcsinh(data[:, j] / COFACTOR)
        if panel is None:
            panel = [
                {
                    "channel": n,
                    "marker": m if t == "marker" else n,
                    "type": t,
                    "cofactor": f"{COFACTOR:g}" if t == "marker" else "",
                    "clustering": "true" if (t == "marker" and m in LINEAGE) else "false",
                    "state": "true" if (t == "marker" and m not in LINEAGE) else "false",
                }
                for n, m, t in zip(common, mk, types)
            ]
        sample_id = row.get("sample_id") or sanitise(Path(row["file"]).stem)
        # Keywords: keep the original ones (minus segment/channel keywords, which are
        # rewritten), then the channel keywords for the kept channels.
        meta = {}
        for k, v in fd.text.items():
            key = "$" + k.upper() if not k.startswith("$") else k.upper()
            if key in SEGMENT_KEYS or re.match(r"^\$P\d+[A-Z]$", key) or key == "$SPILLOVER":
                continue
            if str(v).strip() == "":
                continue
            meta[k] = str(v).strip()
        for j, c in enumerate(common):
            i_src = names.index(c) + 1
            meta[f"p{j + 1}r"] = fd.text.get(f"p{i_src}r", "262144")
            meta[f"p{j + 1}g"] = fd.text.get(f"p{i_src}g", "1")
        meta["pennsieve_prep"] = f"transform=arcsinh;cofactor={COFACTOR:g};fixture=make_fixture.py"
        flat = array("f")
        flat.frombytes(np.ascontiguousarray(data, dtype="<f4").tobytes())
        with open(args.out / "fcs" / f"{sample_id}.fcs", "wb") as fh:
            flowio.create_fcs(fh, flat, common, opt_channel_names=[m if m else None for m in mk],
                              metadata_dict=meta)
        out = {
            "sample_id": sample_id,
            "file": row["file"],
            "batch": row.get("batch", ""),
            "condition": row.get("condition", "NA") or "NA",
            "n_events": data.shape[0],
        }
        for c in extra_cols:
            out[c] = row.get(c, "")
        out_rows.append(out)
        print(f"{sample_id}: {data.shape[0]} events x {data.shape[1]} channels")

    with open(args.out / "metadata" / "samples.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    with open(args.out / "metadata" / "panel.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(panel[0].keys()))
        w.writeheader()
        w.writerows(panel)
    dropped = sorted({c for h in headers.values() for c in h} - set(common))
    print(f"wrote {len(out_rows)} samples to {args.out}; channels dropped (not in every file): {dropped}")


if __name__ == "__main__":
    os.environ.setdefault("PYTHONHASHSEED", "0")
    main()

#!/usr/bin/env python3
"""
Aggregate v4 + overnight BOOST sweep + reproducibility runs into a single
table, derive an updated rule (v2) that uses BOOST when it helps, and dump
JSON + CSV used by build_report_v2.py.

Inputs (relative to this script's repo root, ../data/):
- v4_aggregated.json           (8 datasets x 10 segmenters x {static,auto})
- dataset_meta.json            (GC + bp per dataset)
- dataset_meta_full.json       (chemistry + organism)
- night_runs.csv               (concatenated CSVs from
                                 ~/auto_segmenter/night_runs/out/*/results.csv)

Outputs:
- v5_aggregated.json   combined table keyed by (ds, seg, mode, boost)
- v5_rule_choices.json per-dataset rule pick + measured numbers
- v5_reproducibility.json variance estimates per (ds, seg, mode)
"""
from __future__ import annotations
import csv, json, math
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent
DATA = HERE.parent / "data"

V4 = json.loads((DATA / "v4_aggregated.json").read_text())
REF_META = json.loads((DATA / "dataset_meta.json").read_text())
META = json.loads((DATA / "dataset_meta_full.json").read_text())

# v4 has implicit boost=1.0 for all auto rows.
combined = defaultdict(list)  # key -> list of run dicts
for r in V4:
    key = (r["ds"], r["seg"], r["mode"], 1.0 if r["mode"] == "auto" else None)
    combined[key].append({
        "wall": r["wall"], "ovl01": r["ovl01"], "ovl05": r["ovl05"],
        "d10k": r["d10k"], "d100k": r["d100k"], "n_q": r["n_q"],
        "source": "v4",
    })

def hsap_rescale(f1: float, S: float) -> float:
    """Rescale F1 from full-truth to subsample-truth using fair_hsap.py formula.
    Needed because RB_hsap_sub bench writes raw F1 (against 12974-read truth)
    but the queried fast5 only contains ~2170 reads. S = 12974 / Z."""
    if f1 <= 0 or f1 >= 1:
        return f1
    rho = f1 / (2 - f1)
    rho_new = rho * S
    return 2 * rho_new / (1 + rho_new)

# First pass: figure out hsap_sub default n_q (scaling reference Z).
HSAP_FULL = 12974
hsap_Z = 2170  # default from fair_hsap.py; will be overridden if we have a default row
night_csv = DATA / "night_runs.csv"
if night_csv.exists():
    for row in csv.DictReader(open(night_csv)):
        if row["ds"] == "RB_hsap_sub" and row["seg"] == "default" and row["mode"] == "static":
            try:
                hsap_Z = int(row["n_q"]) or hsap_Z
            except: pass
HSAP_S = HSAP_FULL / hsap_Z

if night_csv.exists():
    for row in csv.DictReader(open(night_csv)):
        try:
            ds, seg, mode = row["ds"], row["seg"], row["mode"]
            # Normalize: static rows always use boost=None so they merge with v4
            # static rows (which have no boost field).
            if mode == "static":
                boost = None
            else:
                boost = float(row["boost"]) if row.get("boost") else 1.0
            ovl01 = float(row["ovl01"]) if row["ovl01"] else 0
            ovl05 = float(row["ovl05"]) if row["ovl05"] else 0
            d10k  = float(row["d10k"])  if row["d10k"]  else 0
            d100k = float(row["d100k"]) if row["d100k"] else 0
            # Apply fair-rescale for RB_hsap_sub so raw night values are
            # directly comparable to v4_aggregated values (which were
            # already rescaled by fix_hsap_eval.py).
            if ds == "RB_hsap_sub":
                ovl01 = hsap_rescale(ovl01, HSAP_S)
                ovl05 = hsap_rescale(ovl05, HSAP_S)
                d10k  = hsap_rescale(d10k,  HSAP_S)
                d100k = hsap_rescale(d100k, HSAP_S)
            combined[(ds, seg, mode, boost)].append({
                "wall": float(row["wall_sec"]),
                "ovl01": ovl01, "ovl05": ovl05,
                "d10k": d10k, "d100k": d100k,
                "n_q":   int(row["n_q"]) if row["n_q"] else 0,
                "source": "night",
            })
        except Exception as e:
            print(f"skip row {row}: {e}")

# Aggregate (median for stability, also keep min/max for the variance plot)
def median(xs):
    s = sorted(xs); n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])

agg = []
for (ds, seg, mode, boost), runs in combined.items():
    walls = [r["wall"] for r in runs if r["wall"] > 0]
    accs  = [r["ovl01"] for r in runs]
    agg.append({
        "ds": ds, "seg": seg, "mode": mode, "boost": boost,
        "n_runs": len(runs),
        "wall_median": median(walls) if walls else 0,
        "wall_min":    min(walls) if walls else 0,
        "wall_max":    max(walls) if walls else 0,
        "ovl01_median": median(accs) if accs else 0,
        "ovl05_median": median([r["ovl05"] for r in runs]) if runs else 0,
        "d10k_median":  median([r["d10k"]  for r in runs]) if runs else 0,
        "d100k_median": median([r["d100k"] for r in runs]) if runs else 0,
        "n_q_median":   median([r["n_q"]   for r in runs]) if runs else 0,
        "wall_cv":      (max(walls)/min(walls) - 1.0) if (walls and min(walls) > 0 and len(walls) > 1) else 0.0,
        "sources": [r["source"] for r in runs],
    })
(DATA / "v5_aggregated.json").write_text(json.dumps(agg, indent=2))


# -------- rule v2 (BOOST-aware) ---------------------------------------------
def best_for(ds: str) -> dict:
    """Return the (seg, mode, boost) row maximizing ovl01_median for this ds.
    Tie-break: lower wall."""
    cands = [r for r in agg if r["ds"] == ds]
    if not cands:
        return {}
    return max(cands, key=lambda r: (r["ovl01_median"], -r["wall_median"]))


def pick_v2(chem: str, ref_mb: float, gc: float | None = None):
    """Rule v2 — data-driven from v4 + overnight BOOST sweep.

    Empirical finding from the overnight RB_ecoli sweep:
      BOOST=1.2 and 1.4 strictly *hurt* accuracy on R10.4.1 top-K segmenters
      (gradient: 0.27 -> 0.08, mad: 0.42 -> 0.10, window: 0.40 -> 0.27),
      contradicting the source comment in revent.c. So we keep BOOST=1.0
      for R10.4.1 and only raise BOOST for plain R10.4 if the D6 sweep
      shows a benefit (decided below from `best_boost_for_d6`)."""
    if chem == "R9.4":
        if ref_mb < 0.1:
            return ("default", "static", None)
        return ("hmm", "static", None)
    if chem == "R10.4":
        # D6 (5.3 Mb): binseg+auto wins.  D7 (2.8 Mb, low GC): default wins,
        # all alternative segmenters land at ovl01 < 0.1.  Split at ref_mb=4.
        if ref_mb < 4:
            return ("default", "static", None)
        return ("binseg", "auto", 1.0)
    if chem == "R10.4.1":
        if ref_mb < 10:
            return ("hmm", "static", None)
        if ref_mb < 1000:
            return ("window", "auto", 1.0)
        return ("default", "static", None)
    return ("default", "static", None)


def lookup(ds, seg, mode, boost):
    """Find aggregate row; fall back to closest boost if exact missing."""
    if mode == "static":
        boost = None
    cands = [r for r in agg if r["ds"] == ds and r["seg"] == seg and r["mode"] == mode]
    if not cands:
        return None
    if boost is None:
        return next((r for r in cands if r["boost"] in (None, 1.0)), cands[0])
    exact = [r for r in cands if r["boost"] == boost]
    if exact:
        return exact[0]
    # nearest boost
    return min(cands, key=lambda r: abs((r["boost"] or 1.0) - boost))


# --- empirically override BOOST for plain R10.4 (D6) if data supports it ---
def best_boost_for(ds: str, seg: str):
    """Among auto rows for (ds, seg), pick the boost with highest median ovl01.
    Tie-break: lower wall."""
    cands = [r for r in agg if r["ds"] == ds and r["seg"] == seg and r["mode"] == "auto"]
    if not cands:
        return None
    return max(cands, key=lambda r: (r["ovl01_median"], -r["wall_median"]))["boost"] or 1.0

BEST_D6_BOOST_BINSEG = best_boost_for("D6", "binseg")
BEST_D6_BOOST_WINDOW = best_boost_for("D6", "window")

rule_rows = []
for ds in sorted({r["ds"] for r in agg}):
    chem = META[ds]["chemistry"]
    ref_mb = META[ds]["ref_mb"]
    gc = REF_META.get(ds, REF_META.get("RB_hsap_full", {})).get("gc", 50.0)
    seg, mode, boost = pick_v2(chem, ref_mb, gc)
    # Data-driven override: for plain R10.4 (D6) use the empirically best BOOST
    if chem == "R10.4" and seg == "binseg" and BEST_D6_BOOST_BINSEG:
        boost = BEST_D6_BOOST_BINSEG
    if chem == "R10.4" and seg == "window" and BEST_D6_BOOST_WINDOW:
        boost = BEST_D6_BOOST_WINDOW
    picked = lookup(ds, seg, mode, boost)
    default = lookup(ds, "default", "static", None)
    oracle = best_for(ds)
    rule_rows.append({
        "ds": ds, "chem": chem, "ref_mb": ref_mb, "gc": gc,
        "rule_seg": seg, "rule_mode": mode, "rule_boost": boost,
        "rule": picked, "default": default, "oracle": oracle,
    })
(DATA / "v5_rule_choices.json").write_text(json.dumps(rule_rows, indent=2, default=str))


# -------- reproducibility / variance summary --------------------------------
repro = []
for r in agg:
    if r["n_runs"] >= 2 and r["wall_min"] > 0:
        repro.append({
            "ds": r["ds"], "seg": r["seg"], "mode": r["mode"], "boost": r["boost"],
            "n_runs": r["n_runs"], "wall_min": r["wall_min"], "wall_med": r["wall_median"],
            "wall_max": r["wall_max"], "wall_pct_spread": r["wall_cv"] * 100,
        })
(DATA / "v5_reproducibility.json").write_text(json.dumps(repro, indent=2))


# -------- console summary ---------------------------------------------------
print(f"aggregated rows: {len(agg)}")
print(f"datasets covered: {sorted({r['ds'] for r in agg})}")
print(f"rule rows: {len(rule_rows)}")
print(f"reproducibility entries (n>=2 runs): {len(repro)}")
print()
print("=== rule decisions vs default & oracle ===")
print(f"{'DS':<12} {'chem':<8} {'ref_mb':>7} {'rule':>22} {'def_acc':>8} "
      f"{'rule_acc':>9} {'orcl_acc':>9} {'speedup':>9} {'dacc':>7}")
for r in rule_rows:
    seg, mode, boost = r["rule_seg"], r["rule_mode"], r["rule_boost"]
    boost_s = "" if boost is None else f"@b{boost}"
    pick_s = f"{seg}/{mode}{boost_s}"
    pacc = r["rule"]["ovl01_median"] if r["rule"] else 0
    pwall = r["rule"]["wall_median"] if r["rule"] else 0
    dacc = r["default"]["ovl01_median"] if r["default"] else 0
    dwall = r["default"]["wall_median"] if r["default"] else 0
    oacc = r["oracle"]["ovl01_median"] if r["oracle"] else 0
    sp = dwall / pwall if pwall else 0
    print(f"{r['ds']:<12} {r['chem']:<8} {r['ref_mb']:>7.1f} {pick_s:>22} "
          f"{dacc:>8.4f} {pacc:>9.4f} {oacc:>9.4f} {sp:>8.2f}x {(pacc-dacc)*100:>+6.1f}")

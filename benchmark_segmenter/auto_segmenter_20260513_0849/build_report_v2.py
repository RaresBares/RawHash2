#!/usr/bin/env python3
"""
Final paper-grade PDF (only graphs).

Reads:
- v5_aggregated.json          (per-cell median + min/max wall, ovl01, etc.)
- v5_rule_choices.json        (per-dataset rule pick with measured stats)
- v5_reproducibility.json     (variance per repeated cell)

Pages:
1. Rule v2 decision table (color-coded segmenter per dataset)
2. Accuracy: default vs rule (v2) vs per-dataset oracle
3. Wall-time per dataset (log) with variance whiskers
4. Speed-up vs default (per dataset)
5. Delta accuracy vs default (per dataset)
6. Auto-K BOOST effect: per (ds, seg), accuracy vs BOOST
7. Auto-K BOOST effect: per (ds, seg), wall vs BOOST
8. Reproducibility: per-cell wall spread (max/min - 1) %
9. All-segmenter Pareto: speed-up vs delta F1, marker = mode, color = family
"""
from __future__ import annotations
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

HERE = Path(__file__).parent
DATA = HERE.parent / "data"

AGG = json.loads((DATA / "v5_aggregated.json").read_text())
CHOICES = json.loads((DATA / "v5_rule_choices.json").read_text())
REPRO = json.loads((DATA / "v5_reproducibility.json").read_text())
META = json.loads((DATA / "dataset_meta_full.json").read_text())

# Sort datasets by (chemistry rank, ref_mb)
def _chem_rank(c): return {"R9.4": 0, "R10.4": 1, "R10.4.1": 2}.get(c, 9)
DATASETS = sorted({r["ds"] for r in AGG}, key=lambda d: (_chem_rank(META[d]["chemistry"]), META[d]["ref_mb"]))

def label(ds):
    m = META[ds]
    return f"{ds}\n{m['organism']}\n{m['chemistry']}"

def get_cell(ds, seg, mode, boost=None):
    if mode == "static": boost = None
    cands = [r for r in AGG if r["ds"] == ds and r["seg"] == seg and r["mode"] == mode]
    if not cands: return None
    if boost is None:
        # prefer the implicit "no boost" row (None or 1.0)
        for c in cands:
            if c["boost"] in (None, 1.0):
                return c
        return cands[0]
    exact = [c for c in cands if c["boost"] == boost]
    return exact[0] if exact else min(cands, key=lambda c: abs((c["boost"] or 1.0) - boost))


SEG_COLOR = {
    "default":   "#666666",
    "scrappie":  "#8c564b",
    "pelt":      "#e377c2",
    "pelt_cuda": "#9467bd",
    "binseg":    "#17becf",
    "hmm":       "#1f77b4",
    "bocd":      "#bcbd22",
    "window":    "#2ca02c",
    "mad":       "#ff7f0e",
    "gradient":  "#d62728",
}


def page_rule_decisions(pdf):
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for i, c in enumerate(CHOICES):
        seg, mode, boost = c["rule_seg"], c["rule_mode"], c["rule_boost"]
        color = SEG_COLOR.get(seg, "#cccccc")
        ax.barh(i, 1, color=color, edgecolor="black", linewidth=0.5)
        b_s = "" if boost is None else f" BOOST={boost}"
        ax.text(0.01, i,
                f"{c['ds']}  ({c['chem']}, {c['ref_mb']:.1f} Mb, GC={c['gc']:.0f}%)"
                f"  ->  {seg}/{mode}{b_s}",
                va="center", ha="left", fontsize=8)
    ax.set_yticks([]); ax.set_xticks([]); ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(CHOICES) - 0.5); ax.invert_yaxis()
    ax.set_title("Rule v2 decisions per dataset")
    h = [plt.Rectangle((0,0),1,1,color=c,label=s) for s,c in SEG_COLOR.items() if s in {c["rule_seg"] for c in CHOICES}]
    ax.legend(handles=h, fontsize=8, loc="lower right")
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_accuracy(pdf):
    fig, ax = plt.subplots(figsize=(12, 6))
    bw = 0.27
    x = np.arange(len(DATASETS))
    d_acc, r_acc, o_acc = [], [], []
    for ds in DATASETS:
        c = next(c for c in CHOICES if c["ds"] == ds)
        d_acc.append(c["default"]["ovl01_median"] if c["default"] else 0)
        r_acc.append(c["rule"]["ovl01_median"]    if c["rule"]    else 0)
        o_acc.append(c["oracle"]["ovl01_median"]  if c["oracle"]  else 0)
    ax.bar(x - bw, d_acc, bw, color="#666666", label="default")
    ax.bar(x,      r_acc, bw, color="#1f77b4", label="rule v2 (Auto-K + BOOST)")
    ax.bar(x + bw, o_acc, bw, color="#2ca02c", label="oracle (best measured)")
    ax.set_xticks(x); ax.set_xticklabels([label(d) for d in DATASETS], fontsize=8)
    ax.set_ylabel("F1 (overlap >= 10%)")
    ax.set_title("F1 (overlap >= 10%): default vs rule v2 vs per-dataset oracle")
    ax.set_ylim(0, max(o_acc) * 1.15); ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_wall(pdf):
    fig, ax = plt.subplots(figsize=(12, 6))
    bw = 0.27
    x = np.arange(len(DATASETS))
    series = {"default": [], "rule": [], "oracle": []}
    err_lo = {"default": [], "rule": [], "oracle": []}
    err_hi = {"default": [], "rule": [], "oracle": []}
    for ds in DATASETS:
        c = next(c for c in CHOICES if c["ds"] == ds)
        for key, cell in [("default", c["default"]), ("rule", c["rule"]), ("oracle", c["oracle"])]:
            w = cell["wall_median"] if cell else 0
            wlo = cell["wall_min"]  if cell else 0
            whi = cell["wall_max"]  if cell else 0
            series[key].append(w)
            err_lo[key].append(max(0, w - wlo) if w > 0 else 0)
            err_hi[key].append(max(0, whi - w) if w > 0 else 0)
    cols = {"default": "#666666", "rule": "#1f77b4", "oracle": "#2ca02c"}
    offs = {"default": -bw, "rule": 0, "oracle": bw}
    for k in series:
        ax.bar(x + offs[k], series[k], bw, color=cols[k], label=k,
               yerr=[err_lo[k], err_hi[k]], capsize=3, error_kw=dict(ecolor="black", elinewidth=0.5))
    ax.set_xticks(x); ax.set_xticklabels([label(d) for d in DATASETS], fontsize=8)
    ax.set_yscale("log"); ax.set_ylabel("Wall time (s, log)")
    ax.set_title("Wall time per dataset (median; whiskers = min/max across runs)")
    ax.legend(fontsize=9); ax.grid(axis="y", which="both", alpha=0.3)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_speedup(pdf):
    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(DATASETS))
    sp = []
    for ds in DATASETS:
        c = next(c for c in CHOICES if c["ds"] == ds)
        dwall = c["default"]["wall_median"] if c["default"] else 1
        rwall = c["rule"]["wall_median"]    if c["rule"]    else 1
        sp.append(dwall / rwall if rwall else 0)
    bars = ax.bar(x, sp, color=["#1f77b4" if s >= 1 else "#d62728" for s in sp])
    ax.axhline(1.0, color="k", linestyle="--", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels([label(d) for d in DATASETS], fontsize=8)
    ax.set_ylabel("Speed-up vs default (default_wall / rule_wall)")
    ax.set_title("Per-dataset wall-time speed-up of rule v2 vs default")
    big = [s for s in sp if s > 0]
    geo = math.exp(np.mean([math.log(s) for s in big])) if big else 0
    ax.text(0.98, 0.95, f"geo-mean = {geo:.2f}x", transform=ax.transAxes,
            ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black"))
    for b, v in zip(bars, sp):
        ax.text(b.get_x() + b.get_width()/2, b.get_height(), f"{v:.2f}x",
                ha="center", va="bottom", fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_accgain(pdf):
    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(DATASETS))
    gain = []
    for ds in DATASETS:
        c = next(c for c in CHOICES if c["ds"] == ds)
        d = c["default"]["ovl01_median"] if c["default"] else 0
        r = c["rule"]["ovl01_median"]    if c["rule"]    else 0
        gain.append((r - d) * 100)
    bars = ax.bar(x, gain, color=["#1f77b4" if g >= 0 else "#d62728" for g in gain])
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels([label(d) for d in DATASETS], fontsize=8)
    ax.set_ylabel("delta F1 (overlap >= 10%) in pp")
    ax.set_title("F1 (overlap) delta: rule v2 - default (positive = rule better)")
    for b, v in zip(bars, gain):
        ax.text(b.get_x()+b.get_width()/2, b.get_height(), f"{v:+.1f}",
                ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_f1_d10k(pdf):
    """Same comparison but using F1 (distance-based, mapping <=10 kb).
    This is the metric the original rawhash paper reports as 'F1' for
    selective-sequencing use cases where only approximate region matters."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    bw = 0.27
    x = np.arange(len(DATASETS))
    d_acc, r_acc, o_acc = [], [], []
    for ds in DATASETS:
        c = next(c for c in CHOICES if c["ds"] == ds)
        cells = [r for r in AGG if r["ds"] == ds]
        oracle_d10k = max(cells, key=lambda r: r["d10k_median"])
        d_acc.append(c["default"]["d10k_median"] if c["default"] else 0)
        r_acc.append(c["rule"]["d10k_median"]    if c["rule"]    else 0)
        o_acc.append(oracle_d10k["d10k_median"])
    ax1.bar(x - bw, d_acc, bw, color="#666666", label="default")
    ax1.bar(x,      r_acc, bw, color="#1f77b4", label="rule v2 (ovl01-optimised)")
    ax1.bar(x + bw, o_acc, bw, color="#2ca02c", label="oracle by d10k")
    ax1.set_xticks(x); ax1.set_xticklabels([label(d) for d in DATASETS], fontsize=8)
    ax1.set_ylabel("F1 (mapping distance <= 10 kb)")
    ax1.set_title("F1 (distance <= 10 kb): default vs rule v2 vs per-dataset oracle\n"
                  "(rule is optimised for F1-overlap; on F1-distance default tends to win)")
    ax1.legend(fontsize=9); ax1.grid(axis="y", alpha=0.3)
    ax1.set_ylim(0, max(o_acc) * 1.1)

    # Delta plot on d10k
    gain = [(r_acc[i] - d_acc[i]) * 100 for i in range(len(DATASETS))]
    bars = ax2.bar(x, gain, color=["#1f77b4" if g >= 0 else "#d62728" for g in gain])
    ax2.axhline(0, color="k", linewidth=0.8)
    ax2.set_xticks(x); ax2.set_xticklabels([label(d) for d in DATASETS], fontsize=8)
    ax2.set_ylabel("delta F1 (d10k) in pp")
    ax2.set_title("F1 (distance) delta: rule v2 - default")
    for b, v in zip(bars, gain):
        ax2.text(b.get_x()+b.get_width()/2, b.get_height(), f"{v:+.1f}",
                ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
    ax2.grid(axis="y", alpha=0.3)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_boost_effect(pdf):
    """For each (ds, seg) where we have multiple boosts, plot accuracy/wall vs boost."""
    # Collect cells grouped by (ds, seg, mode='auto')
    cells = [r for r in AGG if r["mode"] == "auto" and r["seg"] in {"window","mad","binseg","gradient"}]
    keys = sorted({(r["ds"], r["seg"]) for r in cells})
    keys = [k for k in keys if len({r["boost"] for r in cells if r["ds"]==k[0] and r["seg"]==k[1]}) >= 2]
    if not keys:
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))
    for ds, seg in keys:
        rows = sorted([r for r in cells if r["ds"]==ds and r["seg"]==seg],
                       key=lambda r: r["boost"] or 1.0)
        bx = [r["boost"] or 1.0 for r in rows]
        ay = [r["ovl01_median"] for r in rows]
        wy = [r["wall_median"]  for r in rows]
        ax1.plot(bx, ay, "-o", label=f"{ds}/{seg}", color=SEG_COLOR.get(seg, "#666666"))
        ax2.plot(bx, wy, "-o", label=f"{ds}/{seg}", color=SEG_COLOR.get(seg, "#666666"))
    ax1.set_xlabel("RH2_TOPK_AUTO_BOOST"); ax1.set_ylabel("ovl01 (median)")
    ax1.set_title("Auto-K BOOST: accuracy effect (top-K family, R10.4)")
    ax1.grid(alpha=0.3); ax1.legend(fontsize=7, ncol=2)
    ax2.set_xlabel("RH2_TOPK_AUTO_BOOST"); ax2.set_ylabel("wall (s, log)")
    ax2.set_yscale("log")
    ax2.set_title("Auto-K BOOST: wall-time effect")
    ax2.grid(alpha=0.3, which="both"); ax2.legend(fontsize=7, ncol=2)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_reproducibility(pdf):
    fig, ax = plt.subplots(figsize=(11, 5.5))
    rows = sorted(REPRO, key=lambda r: -r["wall_pct_spread"])
    if not rows:
        ax.text(0.5, 0.5, "No reproducibility data (no cells with n_runs >= 2)",
                ha="center", va="center")
    else:
        labels = [f'{r["ds"]}/{r["seg"]}/{r["mode"]}/b={r["boost"] or "-"}' for r in rows[:25]]
        vals = [r["wall_pct_spread"] for r in rows[:25]]
        ns   = [r["n_runs"] for r in rows[:25]]
        y = np.arange(len(labels))
        bars = ax.barh(y, vals,
                       color=["#d62728" if v > 20 else "#ff7f0e" if v > 10 else "#1f77b4" for v in vals])
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("wall-time spread  (max/min - 1) %")
        ax.set_title("Run-to-run reproducibility (top 25 by spread)")
        for b, v, n in zip(bars, vals, ns):
            ax.text(b.get_width(), b.get_y() + b.get_height()/2,
                    f"  {v:.1f}%  (n={n})", va="center", fontsize=7)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def page_pareto(pdf):
    fig, ax = plt.subplots(figsize=(11, 7))
    for ds in DATASETS:
        d = get_cell(ds, "default", "static")
        if not d: continue
        d_wall = d["wall_median"]; d_acc = d["ovl01_median"]
        for r in AGG:
            if r["ds"] != ds or r["seg"] == "default": continue
            sp = d_wall / r["wall_median"] if r["wall_median"] else 0
            dacc = (r["ovl01_median"] - d_acc) * 100
            mk = "o" if r["mode"] == "static" else ("^" if (r["boost"] in (None,1.0)) else "D")
            ax.scatter(sp, dacc, c=SEG_COLOR.get(r["seg"], "#000000"),
                       marker=mk, s=45, alpha=0.7, edgecolors="black", linewidths=0.3)
    ax.axhline(0, color="k", linewidth=0.5); ax.axvline(1, color="k", linewidth=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("Speed-up vs default (log)")
    ax.set_ylabel("delta F1 vs default (ovl01 pt)")
    ax.set_title("All-segmenter Pareto over 8 datasets\no = static, ^ = Auto-K (BOOST=1.0), D = Auto-K + BOOST")
    handles = [plt.Line2D([0],[0], marker="o", color="w", markerfacecolor=c,
                          markersize=8, label=s) for s,c in SEG_COLOR.items() if s != "default"]
    ax.legend(handles=handles, fontsize=7, ncol=2, loc="lower right")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def main():
    out = HERE / "auto_segmenter_v5_report.pdf"
    with PdfPages(out) as pdf:
        page_rule_decisions(pdf)
        page_accuracy(pdf)
        page_wall(pdf)
        page_speedup(pdf)
        page_accgain(pdf)
        page_f1_d10k(pdf)
        page_boost_effect(pdf)
        page_reproducibility(pdf)
        page_pareto(pdf)
    # CSV
    csv_path = HERE / "auto_segmenter_v5_summary.csv"
    with open(csv_path, "w") as fh:
        fh.write("ds,chem,ref_mb,gc,"
                 "default_F1_ovl01,default_F1_d10k,default_wall,"
                 "rule_seg,rule_mode,rule_boost,"
                 "rule_F1_ovl01,rule_F1_d10k,rule_wall,"
                 "oracle_by_ovl01_seg,oracle_by_ovl01_mode,oracle_F1_ovl01,oracle_F1_d10k_for_ovl_oracle,oracle_wall,"
                 "speedup,delta_F1_ovl01_pp,delta_F1_d10k_pp,wall_spread_pct\n")
        for c in CHOICES:
            d = c["default"]; p = c["rule"]; o = c["oracle"]
            sp = (d["wall_median"]/p["wall_median"]) if (p and p["wall_median"]) else 0
            da_ovl = (p["ovl01_median"] - d["ovl01_median"])*100 if (p and d) else 0
            da_d10 = (p["d10k_median"]  - d["d10k_median"]) *100 if (p and d) else 0
            spread = (p["wall_max"]/p["wall_min"] - 1)*100 if (p and p["wall_min"]) else 0
            fh.write(f'{c["ds"]},{c["chem"]},{c["ref_mb"]:.3f},{c["gc"]:.2f},'
                     f'{d["ovl01_median"]:.4f},{d["d10k_median"]:.4f},{d["wall_median"]:.2f},'
                     f'{c["rule_seg"]},{c["rule_mode"]},{c["rule_boost"] or 1.0:.1f},'
                     f'{p["ovl01_median"]:.4f},{p["d10k_median"]:.4f},{p["wall_median"]:.2f},'
                     f'{o["seg"]},{o["mode"]},{o["ovl01_median"]:.4f},{o["d10k_median"]:.4f},{o["wall_median"]:.2f},'
                     f'{sp:.3f},{da_ovl:+.2f},{da_d10:+.2f},{spread:.1f}\n')
    print(f"wrote {out}")
    print(f"wrote {csv_path}")


if __name__ == "__main__":
    main()

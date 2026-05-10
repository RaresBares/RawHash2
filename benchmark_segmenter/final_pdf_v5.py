#!/usr/bin/env python3
"""V5 PDF report: RawHash2 segmenter benchmark with Wilcoxon segmenter added.

Changes vs V4 (final_pdf.py):
  - Adds 'wilcoxon' as the 6th segmenter (sliding-window Mann-Whitney U |z|).
  - New page: Wilcoxon design + synth F1 sweep (0.910 +- 0.007).
  - Per-dataset tables include wilcoxon row.
  - New per-dataset speedup vs default column.
  - Conclusions reflect wilcoxon's real-data behaviour.

Set RESULTS_DIR via env or edit the constant below to run anywhere.
"""
import os
import re
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

RESULTS_DIR = os.environ.get(
    "RESULTS_DIR",
    "/mnt/galactica/rsahleanu/seq_benchmark/results/rh2_100reads",
)
OUT_PDF = os.environ.get(
    "OUT_PDF",
    "/mnt/galactica/rsahleanu/seq_benchmark/results/rawhash2_segmenter_v5.pdf",
)

SEGMENTERS = ["default", "hmm", "pelt", "binseg", "scrappie", "wilcoxon"]
SEG_COLORS = {
    "default": "#2196F3",
    "hmm": "#9C27B0",
    "pelt": "#F44336",
    "binseg": "#FF9800",
    "scrappie": "#4CAF50",
    "wilcoxon": "#00BCD4",
}
SEG_LABELS = {
    "default": "Default\n(t-stat)",
    "hmm": "HMM",
    "pelt": "PELT",
    "binseg": "BinSeg",
    "scrappie": "Scrappie",
    "wilcoxon": "Wilcoxon\n(MW-U)",
}
DS_NAMES = {
    "D1": "D1: SARS-CoV-2 (R9.4)",
    "D2": "D2: E. coli (R9.4)",
    "D3": "D3: Yeast (R9.4)",
    "D4": "D4: Green Algae (R9.4)",
    "D6": "D6: E. coli (R10.4)",
}
DS_GENOME = {
    "D1": "30 kb (virus)",
    "D2": "5 Mb (bacterium)",
    "D3": "12 Mb (yeast, 16 chroms)",
    "D4": "111 Mb (algae, ~17 chroms)",
    "D6": "5 Mb (bacterium, R10.4)",
}


def parse_paf(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return {"total_mappings": 0, "unique_reads": 0, "mapq_ge_60": 0,
                "mean_mapq": 0}
    total = 0
    mapqs = []
    ur = set()
    with open(path) as f:
        for line in f:
            c = line.strip().split("\t")
            if len(c) < 12 or c[5] == "*":
                continue
            try:
                mapq = int(c[11])
            except ValueError:
                continue
            total += 1
            ur.add(c[0])
            mapqs.append(mapq)
    if total == 0:
        return {"total_mappings": 0, "unique_reads": 0, "mapq_ge_60": 0,
                "mean_mapq": 0}
    return {
        "total_mappings": total,
        "unique_reads": len(ur),
        "mapq_ge_60": sum(1 for q in mapqs if q >= 60),
        "mean_mapq": round(sum(mapqs) / len(mapqs), 2),
    }


def parse_time(path):
    if not os.path.exists(path):
        return {"wall_clock_sec": 0, "max_rss_mb": 0}
    d = open(path).read()

    def pw(s):
        p = s.split(":")
        if len(p) == 3:
            return float(p[0]) * 3600 + float(p[1]) * 60 + float(p[2])
        if len(p) == 2:
            return float(p[0]) * 60 + float(p[1])
        return float(s)

    w = re.search(r"Elapsed \(wall clock\) time.*: (.+)", d)
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", d)
    return {
        "wall_clock_sec": round(pw(w.group(1)), 2) if w else 0,
        "max_rss_mb": round(int(m.group(1)) / 1024, 1) if m else 0,
    }


# ---- Load data ----
data = {}
for ds in ["D1", "D2", "D3", "D4", "D6"]:
    ds_dir = os.path.join(RESULTS_DIR, ds)
    if not os.path.isdir(ds_dir):
        continue
    dsd = {"segmenters": {}}
    for seg in SEGMENTERS:
        seg_dir = os.path.join(ds_dir, seg)
        if not os.path.isdir(seg_dir):
            continue
        paf = parse_paf(os.path.join(seg_dir, "mappings.paf"))
        mp = parse_time(os.path.join(seg_dir, "map_time.txt"))
        ip = parse_time(os.path.join(seg_dir, "index_time.txt"))
        acc = {}
        acc_path = os.path.join(seg_dir, "accuracy.json")
        if os.path.exists(acc_path):
            try:
                acc = json.load(open(acc_path))
            except Exception:
                pass
        dsd["segmenters"][seg] = {"paf": paf, "time": mp,
                                  "index_time": ip, "acc": acc}
    data[ds] = dsd


def speedup_vs_default(ds_data, seg):
    """time_default / time_seg. >1 = seg is faster, <1 = seg is slower."""
    if "default" not in ds_data["segmenters"]:
        return None
    if seg not in ds_data["segmenters"]:
        return None
    t_def = ds_data["segmenters"]["default"]["time"]["wall_clock_sec"]
    t_seg = ds_data["segmenters"][seg]["time"]["wall_clock_sec"]
    if t_def <= 0 or t_seg <= 0:
        return None
    return t_def / t_seg


# ---- PDF ----
with PdfPages(OUT_PDF) as pdf:
    # ============ Page 1: Title ============
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    plt.text(0.5, 0.74, "RawHash2 Segmenter\nBenchmark — V5",
             fontsize=36, fontweight="bold", ha="center", va="center",
             transform=fig.transFigure)
    plt.text(0.5, 0.56,
             "Six event segmenters benchmarked end-to-end inside\n"
             "the RawHash2 raw-signal nanopore mapper.",
             fontsize=14, ha="center", va="center",
             transform=fig.transFigure, color="#444")
    plt.text(0.5, 0.42,
             "NEW IN V5: Wilcoxon (Mann-Whitney U sliding-window) segmenter\n"
             "Defaults (w1=5, w2=6, t1=2.5, t2=1.5, ph=0.2) derived from\n"
             "an F1-sweep on synthetic piecewise-constant signals (5 trials,\n"
             "10k samples per trial) — synth F1 = 0.910 +- 0.007.",
             fontsize=11, ha="center", va="center",
             transform=fig.transFigure, color="#555")
    plt.text(0.5, 0.23,
             "Default (t-stat) · HMM · PELT · BinSeg · Scrappie · Wilcoxon\n\n"
             "D1 SARS-CoV-2 · D2 E. coli R9.4 · D3 Yeast · D4 Green Algae · D6 E. coli R10.4",
             fontsize=11, ha="center", va="center",
             transform=fig.transFigure, color="#555")
    plt.text(0.5, 0.06,
             "P&S Arch4Health — Understanding Noise in Genome Sequencers\n"
             "SAFARI Research Group — ETH Zurich — May 2026",
             fontsize=10, ha="center", va="center",
             transform=fig.transFigure, color="#888", style="italic")
    plt.axis("off")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close()

    # ============ Page 2: Background (extended) ============
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.set_title("Background: What is being benchmarked?",
                 fontsize=18, fontweight="bold", loc="left", pad=15)
    text = (
        "\n"
        "Nanopore sequencing measures electrical current as DNA passes through a pore.\n"
        "The raw signal is a noisy step-function: each k-mer in the pore produces a\n"
        "characteristic current level (an 'event').\n"
        "\n"
        "RawHash2 maps raw signals to a reference genome WITHOUT first basecalling\n"
        "(= no base-by-base translation to ACGT text). To make this possible, it must\n"
        "convert the signal into a compact event sequence that can be hashed.\n"
        "\n"
        "THE SEGMENTATION STEP is this conversion: raw samples -> list of events.\n"
        "Each segmenter defines event boundaries differently. If the segmentation is\n"
        "bad, the downstream hash lookups fail and mappings are lost.\n"
        "\n"
        "We ask: how much does the choice of segmenter affect real-world mapping\n"
        "quality, speed, and memory?\n"
        "\n"
        "THE SIX SEGMENTERS IN RAWHASH2:\n"
        "  * default     -- ONT-style dual-window t-statistic (historical default)\n"
        "  * pelt        -- Pruned Exact Linear Time changepoint detection\n"
        "  * binseg      -- Binary (recursive) changepoint detection\n"
        "  * hmm         -- 4-state Hidden Markov Model\n"
        "  * scrappie    -- Scrappie-compatible t-stat with tighter windows\n"
        "  * wilcoxon    -- NEW. Sliding-window Mann-Whitney U |z|-score, replaces\n"
        "                   the Welch t-stat in the same gen_peaks pipeline. Rank-\n"
        "                   based -> robust to heavy-tailed noise / outliers.\n"
        "\n"
        "All six are built into RawHash2 (src/revent.c) and selectable via\n"
        "--segmenter <name>. The downstream signal-to-hash pipeline is identical.\n"
    )
    ax.text(0.02, 0.96, text, fontsize=10.5, family="monospace",
            verticalalignment="top", transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close()

    # ============ Page 3: Wilcoxon design ============
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.set_title("New segmenter: sliding-window Wilcoxon (Mann-Whitney U)",
                 fontsize=16, fontweight="bold", loc="left", pad=15)
    text = (
        "\n"
        "MOTIVATION\n"
        "  The default Welch t-statistic assumes Gaussian per-window data and is\n"
        "  sensitive to outliers (a single spike can shift the mean). Nanopore raw\n"
        "  signal is far from Gaussian under skips, stalls, and pore noise. A rank-\n"
        "  based two-sample test (Mann-Whitney U) trades a small constant-factor\n"
        "  cost for distribution-free robustness.\n"
        "\n"
        "ALGORITHM (revent.c::detect_events_wilcoxon, comp_wilcoxon_z)\n"
        "  At each position i, compare LEFT = signal[i-w .. i] vs RIGHT = signal[i .. i+w]:\n"
        "    U = sum over (l in LEFT, r in RIGHT) of\n"
        "          1.0 if l > r else 0.5 if l == r else 0.0\n"
        "    z = (U - mu) / sigma         mu = w*w / 2,  sigma^2 = w*w(2w+1)/12\n"
        "  Pairwise O(w^2) per position, branch-light, no allocations.\n"
        "  Tie correction omitted (rare for MAD-normalized floats; correction would\n"
        "  shrink sigma^2 -> slightly higher |z|, slightly more peaks).\n"
        "\n"
        "  Then the same dual-window peak-detection pipeline as 'default' (revent.c::\n"
        "  gen_peaks): one short detector + one long detector with threshold-crossing\n"
        "  peak logic, mutual masking. Outputs feed gen_events for per-segment means.\n"
        "\n"
        "F1 SWEEP (synthetic step+Gaussian-noise, mean event length 9, 10000 samples)\n"
        "  Grid: w1 in {2,3,4,5}, w2 in {6,8,10,12,15} (w1<w2), t1 in {1.5..3.5},\n"
        "        t2 in {1.5..3.0}, peak_height in {0.1..0.5}\n"
        "  Boundary tolerance for F1: +-5 samples\n"
        "  5 trials per config, paired (same synth signals across configs)\n"
        "\n"
        "  Optimal mean F1 = 0.910 +- 0.007  @  w1=5  w2=6  t1=2.5  t2=1.5  ph=0.2\n"
        "  Wall time per 10k samples (Python prototype, scipy.stats.mannwhitneyu): 23 ms\n"
        "  C version uses pairwise O(w^2); 5..50x faster than Python for typical w.\n"
        "\n"
        "SOURCES\n"
        "  ~/rawhash2/wilcox_seg/segmenters.py        Python reference\n"
        "  ~/rawhash2/wilcox_seg/sweep_f1.py          synth F1 sweep + diagnostics\n"
        "  ~/rawhash2/src/revent.c                    C port (detect_events_wilcoxon)\n"
        "  ~/rawhash2/src/main.cpp                    --segmenter wilcoxon parsing\n"
        "  ~/rawhash2/src/roptions.h                  RI_SEGMENTER_WILCOXON = 17\n"
    )
    ax.text(0.02, 0.96, text, fontsize=9.5, family="monospace",
            verticalalignment="top", transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close()

    # ============ Page 4: Methodology ============
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.set_title("Methodology", fontsize=18, fontweight="bold", loc="left",
                 pad=15)
    text = (
        "\n"
        "BENCHMARK PROTOCOL (per dataset, per segmenter):\n"
        "  1. Extract 100 reads from FAST5 files (h5py).\n"
        "  2. Build RawHash2 index on reference genome.\n"
        "  3. Map 100 reads -> produce PAF (Pairwise mApping Format).\n"
        "  4. /usr/bin/time -v to record wall-clock and peak RSS.\n"
        "  5. Run minimap2 (-ax map-ont) on basecalled reads -> ground truth PAF.\n"
        "  6. Classify each read via (read_id, reference_id) pair logic.\n"
        "\n"
        "GROUND TRUTH CLASSIFICATION (pair = (read_id, reference_id)):\n"
        "  TP: pair in both rawhash PAF and minimap2 PAF\n"
        "  FP: pair only in rawhash PAF (rawhash claims mapping that truth rejects)\n"
        "  FN: pair only in minimap2 PAF (rawhash missed what truth found)\n"
        "  TN: pair in neither\n"
        "\n"
        "METRICS:\n"
        "  Precision = TP / (TP + FP)        Recall = TP / (TP + FN)\n"
        "  F1        = 2*TP / (2*TP + FP + FN)\n"
        "  Speedup   = wall_clock_default / wall_clock_seg   (>1 = faster than default)\n"
        "\n"
        "HARDWARE: SLURM cluster at ETH Zurich (safari-nexus1), 8 threads, 32 GB RAM.\n"
        "RawHash2 commit: revent.c modified May 2026 to add detect_events_wilcoxon.\n"
        "\n"
        "NOTES:\n"
        "  * R10.4 datasets (D6) require VBZ HDF5 plugin (HDF5_PLUGIN_PATH env).\n"
        "  * D2/D6 reads.fasta not yet available -> no F1 (only mapping stats).\n"
        "  * D3 first 100 reads include MUX-scan calibration reads -> low F1 for all\n"
        "    segmenters; not a deficit specific to any one segmenter.\n"
    )
    ax.text(0.02, 0.96, text, fontsize=10.0, family="monospace",
            verticalalignment="top", transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close()

    # ============ Pages 5-9: Per-dataset tables (with wilcoxon row + speedup) ============
    for ds_key in ["D1", "D2", "D3", "D4", "D6"]:
        if ds_key not in data:
            continue
        ds = data[ds_key]
        if not ds["segmenters"]:
            continue

        fig = plt.figure(figsize=(11, 8.5))
        fig.patch.set_facecolor("white")
        ax = fig.add_subplot(111)
        ax.axis("off")
        title = f"{DS_NAMES[ds_key]}  --  genome: {DS_GENOME[ds_key]}"
        ax.set_title(title, fontsize=14, fontweight="bold", loc="left", pad=10)

        seg_names = [s for s in SEGMENTERS if s in ds["segmenters"]]
        col_labels = ["Segmenter", "Reads\nmapped", "Mean\nMAPQ", "MAPQ>=60",
                      "Precision", "Recall", "F1", "Wall\n(s)", "Speedup\nvs default", "RSS\n(MB)"]
        rows = []
        has_gt = any("precision" in ds["segmenters"][s].get("acc", {})
                     for s in seg_names)
        for s in seg_names:
            sd = ds["segmenters"][s]
            p = sd["paf"]
            t = sd["time"]
            a = sd["acc"]
            sp = speedup_vs_default(ds, s)
            rows.append([
                s.upper(),
                f"{p['unique_reads']}/100",
                f"{p['mean_mapq']:.1f}",
                str(p["mapq_ge_60"]),
                f"{a.get('precision', 0):.3f}" if a else "--",
                f"{a.get('recall', 0):.3f}" if a else "--",
                f"{a.get('f1', 0):.3f}" if a else "--",
                f"{t['wall_clock_sec']:.1f}",
                f"{sp:.2f}x" if sp else "--",
                f"{t['max_rss_mb']:.0f}",
            ])

        table = ax.table(cellText=rows, colLabels=col_labels, loc="upper left",
                         cellLoc="center", colLoc="center",
                         bbox=[0.0, 0.50, 1.0, 0.40])
        table.auto_set_font_size(False)
        table.set_fontsize(9.5)
        for j in range(len(col_labels)):
            table[0, j].set_facecolor("#E3F2FD")
            table[0, j].set_text_props(fontweight="bold")
        for i, s in enumerate(seg_names):
            table[i + 1, 0].set_facecolor(SEG_COLORS[s] + "40")
        # Highlight best F1
        if has_gt:
            f1_vals = [ds["segmenters"][s]["acc"].get("f1", 0) for s in seg_names]
            if any(f1_vals):
                best = f1_vals.index(max(f1_vals))
                table[best + 1, 6].set_facecolor("#C8E6C9")
        # Highlight best wall time
        wt_vals = [ds["segmenters"][s]["time"]["wall_clock_sec"]
                   for s in seg_names]
        wt_nz = [(i, w) for i, w in enumerate(wt_vals) if w > 0]
        if wt_nz:
            best = min(wt_nz, key=lambda x: x[1])[0]
            table[best + 1, 7].set_facecolor("#C8E6C9")
        # Highlight best speedup (excluding default itself)
        sp_vals = [(i, speedup_vs_default(ds, s)) for i, s in enumerate(seg_names)
                   if s != "default"]
        sp_vals = [(i, v) for i, v in sp_vals if v is not None]
        if sp_vals:
            best = max(sp_vals, key=lambda x: x[1])[0]
            table[best + 1, 8].set_facecolor("#C8E6C9")

        # Wilcoxon-specific note per dataset
        wnote = ""
        if "wilcoxon" in ds["segmenters"]:
            w = ds["segmenters"]["wilcoxon"]
            wsp = speedup_vs_default(ds, "wilcoxon")
            wf1 = w["acc"].get("f1") if w["acc"] else None
            if wf1 is not None:
                df1 = ds["segmenters"]["default"]["acc"].get("f1", 0) \
                    if "default" in ds["segmenters"] else 0
                delta = wf1 - df1
                wnote = (
                    f"Wilcoxon on {ds_key}:  F1 = {wf1:.3f}  "
                    f"(default = {df1:.3f}, delta = {delta:+.3f}),  "
                    f"speedup vs default = {wsp:.2f}x"
                    if wsp else f"Wilcoxon on {ds_key}: F1 = {wf1:.3f}"
                )
            else:
                wnote = (f"Wilcoxon on {ds_key}: no ground truth available "
                         f"(reads.fasta missing). Speedup vs default = "
                         f"{wsp:.2f}x" if wsp else f"Wilcoxon on {ds_key}: no GT")
        ax.text(0.02, 0.45, wnote, fontsize=10.5, fontweight="bold",
                color="#00838F",
                verticalalignment="top", transform=ax.transAxes)

        pdf.savefig(fig, bbox_inches="tight")
        plt.close()

    # ============ Page: Wilcoxon F1 per dataset (focused chart) ============
    gt_datasets = [k for k in ["D1", "D3", "D4"]
                   if k in data and "wilcoxon" in data[k]["segmenters"]
                   and "f1" in data[k]["segmenters"]["wilcoxon"].get("acc", {})]
    if gt_datasets:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle("Wilcoxon performance across datasets",
                     fontsize=16, fontweight="bold")

        # Left: F1 wilcoxon vs default per dataset
        ax = axes[0]
        x = np.arange(len(gt_datasets))
        width = 0.35
        f1_w = [data[k]["segmenters"]["wilcoxon"]["acc"].get("f1", 0)
                for k in gt_datasets]
        f1_d = [data[k]["segmenters"]["default"]["acc"].get("f1", 0)
                if "default" in data[k]["segmenters"] else 0
                for k in gt_datasets]
        ax.bar(x - width / 2, f1_d, width, label="default (t-stat)",
               color=SEG_COLORS["default"])
        ax.bar(x + width / 2, f1_w, width, label="wilcoxon (MW-U)",
               color=SEG_COLORS["wilcoxon"])
        ax.set_xticks(x)
        ax.set_xticklabels([DS_NAMES[k].split(":")[0] for k in gt_datasets])
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("F1 score")
        ax.set_title("F1: Wilcoxon vs default")
        ax.legend()
        ax.grid(axis="y", alpha=0.3)
        for i, v in enumerate(f1_d):
            ax.text(x[i] - width / 2, v, f"{v:.2f}",
                    ha="center", va="bottom", fontsize=8)
        for i, v in enumerate(f1_w):
            ax.text(x[i] + width / 2, v, f"{v:.2f}",
                    ha="center", va="bottom", fontsize=8)

        # Right: speedup of wilcoxon vs default per dataset (all 5 if timing
        # available, including no-GT ones)
        ax = axes[1]
        ds_keys_with_w = [k for k in ["D1", "D2", "D3", "D4", "D6"]
                          if k in data and "wilcoxon" in data[k]["segmenters"]]
        sp = [speedup_vs_default(data[k], "wilcoxon") for k in ds_keys_with_w]
        valid = [(k, v) for k, v in zip(ds_keys_with_w, sp) if v is not None]
        if valid:
            ks, vs = zip(*valid)
            colors = ["#43A047" if v >= 1 else "#E53935" for v in vs]
            ax.bar(range(len(ks)), vs, color=colors, alpha=0.85)
            ax.axhline(1.0, color="black", linewidth=1, linestyle="--",
                       label="parity with default")
            ax.set_xticks(range(len(ks)))
            ax.set_xticklabels([DS_NAMES[k].split(":")[0] for k in ks])
            ax.set_ylabel("Speedup vs default  (>1 = faster)")
            ax.set_title("Wilcoxon wall-clock speedup")
            ax.legend()
            ax.grid(axis="y", alpha=0.3)
            for i, v in enumerate(vs):
                ax.text(i, v, f"{v:.2f}x", ha="center", va="bottom",
                        fontsize=9, fontweight="bold")

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close()

    # ============ Page: All-segmenter F1 comparison (datasets with GT) ============
    gt_all = [k for k in ["D1", "D3", "D4"] if k in data and any(
        "precision" in data[k]["segmenters"][s].get("acc", {})
        for s in data[k]["segmenters"])]
    if gt_all:
        n = len(gt_all)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 6), squeeze=False)
        fig.suptitle("F1 across all six segmenters (datasets with ground truth)",
                     fontsize=16, fontweight="bold")
        for idx, ds_key in enumerate(gt_all):
            ax = axes[0, idx]
            ds = data[ds_key]
            seg_names = [s for s in SEGMENTERS if s in ds["segmenters"]
                         and "precision" in ds["segmenters"][s].get("acc", {})]
            x = np.arange(len(seg_names))
            f1s = [ds["segmenters"][s]["acc"].get("f1", 0) for s in seg_names]
            colors = [SEG_COLORS[s] for s in seg_names]
            ax.bar(x, f1s, color=colors, edgecolor="black", linewidth=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels([SEG_LABELS[s] for s in seg_names], fontsize=9)
            ax.set_ylim(0, 1.05)
            ax.set_ylabel("F1 score")
            ax.set_title(DS_NAMES[ds_key])
            ax.grid(axis="y", alpha=0.3)
            for i, v in enumerate(f1s):
                ax.text(x[i], v, f"{v:.3f}", ha="center", va="bottom",
                        fontsize=8, fontweight="bold")
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close()

    # ============ Page: Speed vs accuracy scatter ============
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Speed vs Accuracy", fontsize=16, fontweight="bold")

    ax = axes[0]
    for ds_key in ["D1", "D2", "D3", "D4", "D6"]:
        if ds_key not in data:
            continue
        for seg in data[ds_key]["segmenters"]:
            p = data[ds_key]["segmenters"][seg]["paf"]
            t = data[ds_key]["segmenters"][seg]["time"]
            if t["wall_clock_sec"] > 0:
                ax.scatter(t["wall_clock_sec"], p["unique_reads"],
                           s=100, color=SEG_COLORS[seg], edgecolor="black",
                           alpha=0.7, linewidth=0.5)
    for seg in SEGMENTERS:
        ax.scatter([], [], s=100, color=SEG_COLORS[seg], edgecolor="black",
                   label=seg)
    ax.set_xscale("log")
    ax.set_xlabel("Wall clock time (s, log scale)")
    ax.set_ylabel("Unique reads mapped (out of 100)")
    ax.set_title("All datasets: reads mapped vs time")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    if gt_all:
        markers = {"D1": "o", "D3": "s", "D4": "^"}
        for ds_key in gt_all:
            ds = data[ds_key]
            for seg in ds["segmenters"]:
                a = ds["segmenters"][seg]["acc"]
                t = ds["segmenters"][seg]["time"]
                if "f1" in a and t["wall_clock_sec"] > 0:
                    ax.scatter(t["wall_clock_sec"], a["f1"],
                               s=120, color=SEG_COLORS[seg], edgecolor="black",
                               linewidth=0.8, marker=markers.get(ds_key, "o"),
                               alpha=0.85)
                    ax.annotate(seg[:3], (t["wall_clock_sec"], a["f1"]),
                                textcoords="offset points", xytext=(7, 3),
                                fontsize=8)
        for ds_key, m in markers.items():
            if ds_key in gt_all:
                ax.scatter([], [], marker=m, s=120, c="gray",
                           edgecolor="black",
                           label=DS_NAMES[ds_key].split(":")[0])
    ax.set_xscale("log")
    ax.set_xlabel("Wall clock time (s, log scale)")
    ax.set_ylabel("F1 score")
    ax.set_title("F1 vs speed (datasets with ground truth)")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close()

    # ============ Page: Conclusions ============
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.set_title("Conclusions & V5 Findings",
                 fontsize=18, fontweight="bold", loc="left", pad=15)

    # Compute summary numbers
    wf1_lines = []
    wsp_lines = []
    for ds_key in ["D1", "D2", "D3", "D4", "D6"]:
        if ds_key not in data:
            continue
        if "wilcoxon" not in data[ds_key]["segmenters"]:
            continue
        w = data[ds_key]["segmenters"]["wilcoxon"]
        sp = speedup_vs_default(data[ds_key], "wilcoxon")
        if "f1" in w["acc"]:
            df1 = data[ds_key]["segmenters"].get("default", {}).get("acc", {}).get("f1", 0)
            wf1_lines.append(
                f"  {ds_key:3s}  F1 = {w['acc']['f1']:.3f}  "
                f"(default = {df1:.3f}, delta = {w['acc']['f1']-df1:+.3f})"
            )
        if sp is not None:
            wsp_lines.append(f"  {ds_key:3s}  speedup = {sp:.2f}x")

    text = (
        "\n"
        "WILCOXON SEGMENTER -- DATASET-LEVEL RESULTS\n"
        "===========================================================================\n"
        "\n"
        "Per-dataset F1 (where ground truth available):\n"
        + ("\n".join(wf1_lines) if wf1_lines else "  (no F1 available)") +
        "\n\n"
        "Per-dataset speedup vs default (>1 = wilcoxon faster):\n"
        + ("\n".join(wsp_lines) if wsp_lines else "  (no timing available)") +
        "\n\n"
        "KEY FINDINGS\n"
        "===========================================================================\n"
        "\n"
        "1. Wilcoxon is competitive on hard genomes, weaker on easy ones.\n"
        "   On D1 (SARS-CoV-2, 30 kb): F1 = 0.816 vs default 0.953 -- worse.\n"
        "   On D3 (yeast, 12 Mb):       F1 = 0.417 vs default 0.957 -- much worse.\n"
        "   On D4 (algae, 111 Mb):      F1 = 0.344 vs default 0.265 -- WINS by +0.08\n"
        "   D4 is the hardest dataset (most chromosomes, lowest absolute F1 across\n"
        "   all segmenters); the rank-based statistic seems to handle the noise\n"
        "   floor / outlier regime there better than the t-statistic.\n"
        "\n"
        "2. Wilcoxon is slower than default on every dataset (speedup < 1x).\n"
        "   Range: 0.25x (D1) to 0.87x (D2). Pairwise O(n*w^2) inner loop costs\n"
        "   ~2-4x of the t-stat on the segmentation step, but mapping wall-time\n"
        "   is dominated by indexing + chaining + I/O so the end-to-end gap stays\n"
        "   in the 1.1-4x range, never above parity.\n"
        "\n"
        "3. Default/PELT win on D1, D2, D3; Wilcoxon wins on D4.\n"
        "   No segmenter dominates everywhere. Choice depends on signal regime:\n"
        "   stationary Gaussian -> t-stat; heavy-tailed / drifty -> Wilcoxon.\n"
        "\n"
        "4. BinSeg still broken (F1 ~ 0.02-0.06) -- unchanged from V4.\n"
        "\n"
        "FUTURE WORK\n"
        "===========================================================================\n"
        "  * Per-dataset Wilcoxon hyperparameter sweep on real signals (the synth-\n"
        "    derived defaults are not optimal for real nanopore data).\n"
        "  * Try a robust-scaler-normalized Wilcoxon (median+MAD per window) to\n"
        "    eliminate slow drift before ranking.\n"
        "  * Complete D2/D6 ground truth once SRA reads.fasta retrievals finish.\n"
        "  * Compare event count distributions: does Wilcoxon under-segment in\n"
        "    no-event regions and over-segment in transitions?\n"
    )
    ax.text(0.02, 0.96, text, fontsize=9.5, family="monospace",
            verticalalignment="top", transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close()

print(f"PDF saved: {OUT_PDF}")

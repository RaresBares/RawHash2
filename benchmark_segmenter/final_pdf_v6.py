#!/usr/bin/env python3
"""V6 PDF report: RawHash2 segmenter benchmark.

V6 changes vs V5:
  - Wilcoxon now uses DEFAULT-ALIGNED config: same windows as t-stat default
    (w1=3, w2=9, ph=0.4) with thresholds rescaled (t1=1.8, t2=3.0) to
    Wilcoxon's bounded |z| range.
  - Per-dataset F1 measurements with this corrected config.
  - New page: why direct 1:1 transfer from t-stat fails (bounded |z|).
  - Within-plateau noise sigma^2 measured on real D1 reads.

Inputs (env vars):
  RESULTS_DIR    base bench dir
  REBENCH_DIR    rebench dir for default-aligned wilcoxon
  OUT_PDF        output path
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

RESULTS_DIR = os.environ.get("RESULTS_DIR",
    "/mnt/galactica/rsahleanu/seq_benchmark/results/rh2_100reads")
REBENCH_DIR = os.environ.get("REBENCH_DIR",
    "/mnt/galactica/rsahleanu/seq_benchmark/results/wilcox_dflt_aligned")
OUT_PDF = os.environ.get("OUT_PDF",
    "/mnt/galactica/rsahleanu/seq_benchmark/results/rawhash2_segmenter_v6.pdf")

SEGMENTERS = ["default", "hmm", "pelt", "binseg", "scrappie", "wilcoxon"]
SEG_COLORS = {"default": "#2196F3", "hmm": "#9C27B0", "pelt": "#F44336",
              "binseg": "#FF9800", "scrappie": "#4CAF50", "wilcoxon": "#00BCD4"}
SEG_LABELS = {"default": "Default\n(t-stat)", "hmm": "HMM", "pelt": "PELT",
              "binseg": "BinSeg", "scrappie": "Scrappie",
              "wilcoxon": "Wilcoxon\n(default-aligned)"}
DS_NAMES = {"D1": "D1: SARS-CoV-2 (R9.4)", "D2": "D2: E. coli (R9.4)",
            "D3": "D3: Yeast (R9.4)", "D4": "D4: Green Algae (R9.4)",
            "D6": "D6: E. coli (R10.4)"}
DS_GENOME = {"D1": "30 kb (virus)", "D2": "5 Mb (bacterium)",
             "D3": "12 Mb (yeast, 16 chroms)",
             "D4": "111 Mb (algae, ~17 chroms)",
             "D6": "5 Mb (bacterium, R10.4)"}


def parse_paf(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return {"total_mappings": 0, "unique_reads": 0, "mapq_ge_60": 0, "mean_mapq": 0}
    total = 0; mapqs = []; ur = set()
    with open(path) as f:
        for line in f:
            c = line.strip().split("\t")
            if len(c) < 12 or c[5] == "*": continue
            try: mapq = int(c[11])
            except ValueError: continue
            total += 1; ur.add(c[0]); mapqs.append(mapq)
    if total == 0:
        return {"total_mappings": 0, "unique_reads": 0, "mapq_ge_60": 0, "mean_mapq": 0}
    return {"total_mappings": total, "unique_reads": len(ur),
            "mapq_ge_60": sum(1 for q in mapqs if q >= 60),
            "mean_mapq": round(sum(mapqs) / len(mapqs), 2)}


def parse_time(path):
    if not os.path.exists(path):
        return {"wall_clock_sec": 0, "max_rss_mb": 0}
    d = open(path).read()
    def pw(s):
        p = s.split(":")
        if len(p) == 3: return float(p[0])*3600 + float(p[1])*60 + float(p[2])
        if len(p) == 2: return float(p[0])*60 + float(p[1])
        return float(s)
    w = re.search(r"Elapsed \(wall clock\) time.*: (.+)", d)
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", d)
    return {"wall_clock_sec": round(pw(w.group(1)), 2) if w else 0,
            "max_rss_mb": round(int(m.group(1)) / 1024, 1) if m else 0}


# ---- Load bench data ----
data = {}
for ds in ["D1", "D2", "D3", "D4", "D6"]:
    ds_dir = os.path.join(RESULTS_DIR, ds)
    if not os.path.isdir(ds_dir): continue
    dsd = {"segmenters": {}}
    for seg in SEGMENTERS:
        seg_dir = None
        if seg == "wilcoxon":
            cand = os.path.join(REBENCH_DIR, ds, "wilcoxon")
            if os.path.isdir(cand): seg_dir = cand
        if seg_dir is None:
            seg_dir = os.path.join(ds_dir, seg)
        if not os.path.isdir(seg_dir): continue
        paf = parse_paf(os.path.join(seg_dir, "mappings.paf"))
        mp = parse_time(os.path.join(seg_dir, "map_time.txt"))
        ip = parse_time(os.path.join(seg_dir, "index_time.txt"))
        acc = {}
        acc_path = os.path.join(seg_dir, "accuracy.json")
        if os.path.exists(acc_path):
            try: acc = json.load(open(acc_path))
            except: pass
        dsd["segmenters"][seg] = {"paf": paf, "time": mp, "index_time": ip,
                                  "acc": acc, "_dir": seg_dir}
    data[ds] = dsd


def speedup(ds_data, seg):
    if "default" not in ds_data["segmenters"] or seg not in ds_data["segmenters"]:
        return None
    t_def = ds_data["segmenters"]["default"]["time"]["wall_clock_sec"]
    t_seg = ds_data["segmenters"][seg]["time"]["wall_clock_sec"]
    if t_def <= 0 or t_seg <= 0: return None
    return t_def / t_seg


WILCOX_CFG = "w1=3, w2=9, t1=1.8, t2=3.0, ph=0.4"

# ---- PDF ----
with PdfPages(OUT_PDF) as pdf:
    # ============ Page 1: Title ============
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    plt.text(0.5, 0.74, "RawHash2 Segmenter\nBenchmark — V6",
             fontsize=36, fontweight="bold", ha="center", va="center",
             transform=fig.transFigure)
    plt.text(0.5, 0.55,
             "Wilcoxon segmenter realigned to default's window structure.\n"
             "Same w1=3, w2=9, peak_height=0.4 as t-stat default;\n"
             "thresholds rescaled to Wilcoxon's bounded |z|: t1=1.8, t2=3.0.",
             fontsize=13, ha="center", va="center",
             transform=fig.transFigure, color="#444")
    plt.text(0.5, 0.39,
             "RESULT: Wilcoxon recovers most of default's accuracy on D1\n"
             "(0.909 vs 0.953) and beats default on the hard D4 algae\n"
             "genome (0.374 vs 0.265).",
             fontsize=12, ha="center", va="center",
             transform=fig.transFigure, color="#00838F", fontweight="bold")
    plt.text(0.5, 0.22,
             "Default · HMM · PELT · BinSeg · Scrappie · Wilcoxon (default-aligned)\n\n"
             "D1 SARS-CoV-2 · D2 E. coli R9.4 · D3 Yeast · D4 Green Algae · D6 E. coli R10.4",
             fontsize=11, ha="center", va="center",
             transform=fig.transFigure, color="#555")
    plt.text(0.5, 0.06,
             "P&S Arch4Health — Understanding Noise in Genome Sequencers\n"
             "SAFARI Research Group — ETH Zurich — May 2026",
             fontsize=10, ha="center", va="center",
             transform=fig.transFigure, color="#888", style="italic")
    plt.axis("off")
    pdf.savefig(fig, bbox_inches="tight"); plt.close()

    # ============ Page 2: Why default-aligned ============
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111); ax.axis("off")
    ax.set_title("Why default-aligned Wilcoxon (and why direct 1:1 fails)",
                 fontsize=15, fontweight="bold", loc="left", pad=15)
    text = (
        "\n"
        "QUESTION: can we just use the t-stat default's parameters with\n"
        "a Wilcoxon test underneath?\n"
        "\n"
        "  default t-stat config:  w1=3, w2=9, t1=4.0, t2=3.5, ph=0.4\n"
        "\n"
        "ANSWER: not directly, because Wilcoxon's |z| is BOUNDED:\n"
        "\n"
        "  Welch's |t|:    asymptotically unbounded (variance in denominator\n"
        "                  can be arbitrarily small)\n"
        "  Wilcoxon |z|:   bounded by |z|_max ≈ sqrt(3w/2)\n"
        "                  -> w=3:  max |z| = 1.97   (t1=4.0 unreachable!)\n"
        "                  -> w=9:  max |z| = 3.58   (t2=3.5 borderline)\n"
        "                  -> w=15: max |z| = 4.60\n"
        "\n"
        "  Setting t1=4.0 with w=3 means the short detector NEVER triggers.\n"
        "  The asymptotic-normal approximation only kicks in at w >= ~20,\n"
        "  but nanopore demands small windows (~k-mer length / sample rate).\n"
        "\n"
        "RESCALING via equivalent p-value:\n"
        "\n"
        "  default t-stat |t|=4.0 (df=2(w-1)=4)  -> p ≈ 0.016\n"
        "    -> equivalent Wilcoxon |z| (two-sided) ≈ 2.41\n"
        "    -> still > 1.97 = max for w=3.  Best feasible approximation: 1.8.\n"
        "\n"
        "  default t-stat |t|=3.5 (df=16)        -> p ≈ 0.003\n"
        "    -> equivalent Wilcoxon |z| ≈ 3.00\n"
        "    -> achievable with w=9 (max |z| = 3.58).\n"
        "\n"
        "WILCOXON DEFAULT-ALIGNED CONFIG:  " + WILCOX_CFG + "\n"
        "\n"
        "PEAK_HEIGHT (0.4) is dimensionless on |z| just as on |t| -- preserved.\n"
        "Window sizes (w1=3, w2=9) preserved -- same scale as t-stat default.\n"
    )
    ax.text(0.02, 0.96, text, fontsize=10.5, family="monospace",
            verticalalignment="top", transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight"); plt.close()

    # ============ Page 3: Headline F1 chart ============
    gt_keys = [k for k in ["D1", "D3", "D4"] if k in data and
               "wilcoxon" in data[k]["segmenters"] and
               "f1" in data[k]["segmenters"]["wilcoxon"].get("acc", {})]
    if gt_keys:
        fig, ax = plt.subplots(1, 1, figsize=(11, 6))
        x = np.arange(len(gt_keys))
        width = 0.35
        f1_d = [data[k]["segmenters"]["default"]["acc"].get("f1", 0) for k in gt_keys]
        f1_w = [data[k]["segmenters"]["wilcoxon"]["acc"]["f1"] for k in gt_keys]
        ax.bar(x - width/2, f1_d, width, color=SEG_COLORS["default"],
               label="Default (t-stat)")
        ax.bar(x + width/2, f1_w, width, color=SEG_COLORS["wilcoxon"],
               label=f"Wilcoxon ({WILCOX_CFG})")
        ax.set_xticks(x)
        ax.set_xticklabels([DS_NAMES[k] for k in gt_keys])
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("F1 score (mapping vs minimap2 truth)")
        ax.set_title("V6: Wilcoxon (default-aligned) vs Default t-stat",
                     fontsize=14, fontweight="bold")
        ax.legend(fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        for i, v in enumerate(f1_d):
            ax.text(x[i] - width/2, v, f"{v:.3f}", ha="center", va="bottom", fontsize=9)
        for i, v in enumerate(f1_w):
            color = "#00838F" if f1_w[i] >= f1_d[i] else "#999"
            weight = "bold" if f1_w[i] >= f1_d[i] else "normal"
            ax.text(x[i] + width/2, v, f"{v:.3f}", ha="center", va="bottom",
                    fontsize=9, fontweight=weight, color=color)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close()

    # ============ Page 4: Per-dataset tables ============
    for ds_key in ["D1", "D2", "D3", "D4", "D6"]:
        if ds_key not in data: continue
        ds = data[ds_key]
        if not ds["segmenters"]: continue

        fig = plt.figure(figsize=(11, 8.5))
        fig.patch.set_facecolor("white")
        ax = fig.add_subplot(111); ax.axis("off")
        title = f"{DS_NAMES[ds_key]}  --  genome: {DS_GENOME[ds_key]}"
        ax.set_title(title, fontsize=14, fontweight="bold", loc="left", pad=10)

        seg_names = [s for s in SEGMENTERS if s in ds["segmenters"]]
        col_labels = ["Segmenter", "Reads\nmapped", "Mean\nMAPQ", "MAPQ>=60",
                      "Precision", "Recall", "F1", "Wall\n(s)",
                      "Speedup\nvs default", "RSS\n(MB)"]
        rows = []
        has_gt = any("precision" in ds["segmenters"][s].get("acc", {})
                     for s in seg_names)
        for s in seg_names:
            sd = ds["segmenters"][s]
            p = sd["paf"]; t = sd["time"]; a = sd["acc"]
            sp = speedup(ds, s)
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
        table.auto_set_font_size(False); table.set_fontsize(9.5)
        for j in range(len(col_labels)):
            table[0, j].set_facecolor("#E3F2FD")
            table[0, j].set_text_props(fontweight="bold")
        for i, s in enumerate(seg_names):
            table[i + 1, 0].set_facecolor(SEG_COLORS[s] + "40")
        if has_gt:
            f1_vals = [ds["segmenters"][s]["acc"].get("f1", 0) for s in seg_names]
            if any(f1_vals):
                best = f1_vals.index(max(f1_vals))
                table[best + 1, 6].set_facecolor("#C8E6C9")
        wt_vals = [ds["segmenters"][s]["time"]["wall_clock_sec"] for s in seg_names]
        wt_nz = [(i, w) for i, w in enumerate(wt_vals) if w > 0]
        if wt_nz:
            best = min(wt_nz, key=lambda x: x[1])[0]
            table[best + 1, 7].set_facecolor("#C8E6C9")

        if "wilcoxon" in ds["segmenters"]:
            w = ds["segmenters"]["wilcoxon"]
            wsp = speedup(ds, "wilcoxon")
            wf1 = w["acc"].get("f1") if w["acc"] else None
            wnote = f"Wilcoxon (V6 default-aligned, {WILCOX_CFG}):\n"
            if wf1 is not None:
                df1 = ds["segmenters"].get("default", {}).get("acc", {}).get("f1", 0)
                delta = wf1 - df1
                wnote += (f"  F1 = {wf1:.3f}  (default {df1:.3f}, delta = {delta:+.3f})"
                          f"{', WINS' if delta>0 else ''},  "
                          f"speedup = {wsp:.2f}x" if wsp else f"  F1 = {wf1:.3f}")
            elif wsp:
                wnote += f"  No GT for this dataset; speedup vs default = {wsp:.2f}x"
            ax.text(0.02, 0.45, wnote, fontsize=10.5, fontweight="bold",
                    color="#00838F", verticalalignment="top",
                    transform=ax.transAxes)
        pdf.savefig(fig, bbox_inches="tight"); plt.close()

    # ============ Page: Within-plateau noise sigma^2 (CRLB measurement) ============
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111); ax.axis("off")
    ax.set_title("Within-plateau additive-noise sigma^2 (for CRLB)",
                 fontsize=15, fontweight="bold", loc="left", pad=15)
    text = (
        "\n"
        "MODEL\n"
        "  sig[t] = mu_kmer + noise[t],     noise[t] ~ N(0, sigma^2)\n"
        "\n"
        "MEASUREMENT (10 reads from D1 SARS-CoV-2 / R9.4 fast5)\n"
        "\n"
        "  Procedure:\n"
        "    1. Read raw signal, convert DAC -> pA via channel calibration.\n"
        "    2. Run default t-stat segmenter to identify plateaus.\n"
        "    3. Per plateau (>=3 samples), compute residual variance.\n"
        "    4. Aggregate: mean and median across plateaus.\n"
        "\n"
        "  Result (10 reads, ~430 plateaus per read):\n"
        "                          MEAN          MEDIAN (robust)\n"
        "    sigma^2 (raw pA^2):   23.10 pA^2    13.34 pA^2\n"
        "    sigma   (raw pA):      4.81 pA       3.65 pA\n"
        "    sigma^2 (MAD-norm):    0.158         0.083\n"
        "    sigma   (MAD-norm):    0.40          0.29\n"
        "    Channel MAD:                         12.5 pA\n"
        "\n"
        "RECOMMENDATION FOR CRLB\n"
        "  Use sigma^2 ≈ 13 pA^2 (median, robust to pore-stall outliers)\n"
        "  Use sigma   ≈ 3.65 pA\n"
        "\n"
        "  The mean is inflated by occasional reads with stalled pores; one read\n"
        "  in the sample had sigma^2=109 pA^2 vs ~12 pA^2 typical -> the median\n"
        "  is the cleaner estimate of the true sensor noise.\n"
        "\n"
        "  This measurement is on R9.4 chemistry; R10.4 (D6) noise floor differs;\n"
        "  re-measure on D6 fast5 if the CRLB analysis covers R10.4.\n"
    )
    ax.text(0.02, 0.96, text, fontsize=10.5, family="monospace",
            verticalalignment="top", transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight"); plt.close()

    # ============ Page: Conclusions ============
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111); ax.axis("off")
    ax.set_title("V6 Conclusions", fontsize=18, fontweight="bold",
                 loc="left", pad=15)

    summary = []
    for ds_key in ["D1", "D2", "D3", "D4", "D6"]:
        if ds_key not in data or "wilcoxon" not in data[ds_key]["segmenters"]:
            continue
        w = data[ds_key]["segmenters"]["wilcoxon"]
        sp = speedup(data[ds_key], "wilcoxon")
        if "f1" in w["acc"]:
            df1 = data[ds_key]["segmenters"].get("default", {}).get("acc", {}).get("f1", 0)
            wins = " WINS" if w["acc"]["f1"] > df1 else ""
            summary.append(
                f"  {ds_key:3s}  F1 = {w['acc']['f1']:.3f}  "
                f"(default {df1:.3f}, delta = {w['acc']['f1']-df1:+.3f}){wins:5s}  "
                f"speedup = {sp:.2f}x" if sp else
                f"  {ds_key:3s}  F1 = {w['acc']['f1']:.3f}{wins}"
            )
        elif sp:
            summary.append(f"  {ds_key:3s}  no GT,  speedup = {sp:.2f}x")

    text = (
        "\n"
        "WILCOXON (V6 DEFAULT-ALIGNED)  vs  DEFAULT (t-stat)\n"
        "===========================================================================\n"
        "\n"
        + ("\n".join(summary) if summary else "  (no data)") +
        "\n\n"
        "TAKEAWAYS\n"
        "===========================================================================\n"
        "\n"
        "1. Default-aligned Wilcoxon is competitive again.\n"
        "   By keeping w1=3, w2=9, ph=0.4 (same as t-stat default) and only\n"
        "   rescaling thresholds to Wilcoxon's bounded |z| range (1.8, 3.0),\n"
        "   Wilcoxon F1 jumps from synth-baked ~0.81 to 0.91 on D1 -- closing\n"
        "   most of the gap to default's 0.95.\n"
        "\n"
        "2. Wilcoxon WINS on D4 (Green Algae, hard genome).\n"
        "   F1 = 0.374 vs default 0.265 (+0.11 absolute). On the noisiest /\n"
        "   most heterogeneous genome we tested, the rank-based statistic\n"
        "   handles the signal regime better than the t-stat.\n"
        "\n"
        "3. Default still wins on the well-behaved genomes (D1, D3).\n"
        "   On clean Gaussian-like signal the t-stat's parametric leverage\n"
        "   gives it the edge.\n"
        "\n"
        "4. Methodological lesson: do NOT use synth boundary-F1 to pick params.\n"
        "   The first Wilcoxon defaults (V5: w1=5, w2=6, t1=2.5, t2=1.5, ph=0.2)\n"
        "   were synth-derived -- F1=0.910 on synth, only 0.816 on real D1.\n"
        "   Default-aligned config (V6: w1=3, w2=9, t1=1.8, t2=3.0, ph=0.4) is\n"
        "   a +0.09 jump on real data.\n"
        "\n"
        "FUTURE WORK\n"
        "===========================================================================\n"
        "  * Per-dataset fine-tuning around (3, 9, 1.8, 3.0, 0.4).\n"
        "  * Add tie correction to comp_wilcoxon_z (currently dropped).\n"
        "  * Measure within-plateau sigma^2 on D6 (R10.4) for that CRLB.\n"
        "  * D4 ground-truth F1 is uniformly low (default 0.265, wilcoxon 0.374);\n"
        "    investigate whether it's a fundamental signal/index issue rather\n"
        "    than a segmenter issue.\n"
    )
    ax.text(0.02, 0.96, text, fontsize=9.5, family="monospace",
            verticalalignment="top", transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches="tight"); plt.close()

print(f"PDF saved: {OUT_PDF}")

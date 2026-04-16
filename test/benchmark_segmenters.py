#!/usr/bin/env python3
"""
Segmenter Benchmark for RawHash2
Runs RawHash2 with different event segmenters on multiple datasets,
collects mapping metrics and segmentation quality, and generates a PDF report.
"""

import os
import sys
import json
import time
import glob
import subprocess
import numpy as np
from collections import defaultdict

# VBZ plugin for R10.4 FAST5
os.environ["HDF5_PLUGIN_PATH"] = os.path.expanduser(
    "~/.local/lib/python3.10/site-packages/vbz_h5py_plugin/lib")

import h5py

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

# =============================================================================
# Configuration
# =============================================================================

RAWHASH2 = os.path.expanduser("~/rawhash2/bin/rawhash2")
DATA_DIR = os.path.expanduser("~/rawhash2/test/data")
PORE_R94 = os.path.expanduser("~/rawhash2/extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model")
PORE_R104 = os.path.expanduser("~/rawhash2/extern/local_kmer_models/uncalled_r1041_model_only_means.txt")
OUTPUT_DIR = os.path.expanduser("~/rawhash2/test/benchmark_results")
MAX_READS = 50  # Only process first 50 reads

DATASETS = {
    "d1_sars-cov-2_r94": {"name": "D1: SARS-CoV-2 (R9.4)", "pore": PORE_R94, "preset": "viral", "extra": "", "size": "small"},
    "d7_saureus_r104":   {"name": "D7: S. aureus (R10.4)",  "pore": PORE_R104, "preset": "sensitive", "extra": "--r10", "size": "small"},
    "d2_ecoli_r94":      {"name": "D2: E. coli (R9.4)",     "pore": PORE_R94, "preset": "sensitive", "extra": "", "size": "medium"},
    "d6_ecoli_r104":     {"name": "D6: E. coli (R10.4)",    "pore": PORE_R104, "preset": "sensitive", "extra": "--r10", "size": "medium"},
    "d3_yeast_r94":      {"name": "D3: Yeast (R9.4)",       "pore": PORE_R94, "preset": "sensitive", "extra": "", "size": "large"},
    "d4_green_algae_r94":{"name": "D4: Green Algae (R9.4)", "pore": PORE_R94, "preset": "sensitive", "extra": "", "size": "large"},
}

SEGMENTERS = ["default", "hmm", "pelt", "binseg", "window"]
SEGMENTER_LABELS = {
    "default": "ONT t-stat (default)",
    "hmm": "HMM (Viterbi)",
    "pelt": "PELT (changepoint)",
    "binseg": "BinSeg (recursive)",
    "window": "Window (sliding)",
}

SEGMENTER_COLORS = {
    "default": "#1f77b4",
    "hmm": "#ff7f0e",
    "pelt": "#2ca02c",
    "binseg": "#d62728",
    "window": "#9467bd",
}

SEGMENTER_CATEGORIES = {
    "Statistical (t-test)": ["default"],
    "Probabilistic (HMM)": ["hmm"],
    "Changepoint Detection": ["pelt", "binseg"],
    "Window-based": ["window"],
}

# =============================================================================
# Helper functions
# =============================================================================

def prepare_fast5_subset(dataset_id, max_reads=MAX_READS):
    """Create a subset directory with only first max_reads FAST5 files."""
    fast5_dir = os.path.join(DATA_DIR, dataset_id, "fast5_files")
    subset_dir = os.path.join(DATA_DIR, dataset_id, "fast5_subset")

    if not os.path.exists(fast5_dir):
        print(f"  WARNING: {fast5_dir} not found")
        return None

    os.makedirs(subset_dir, exist_ok=True)

    # Check if subset already exists with enough files
    existing = glob.glob(os.path.join(subset_dir, "*.fast5"))
    if len(existing) >= max_reads:
        return subset_dir

    # Clear and rebuild
    for f in existing:
        os.remove(f)

    fast5_files = sorted(glob.glob(os.path.join(fast5_dir, "*.fast5")))[:max_reads]
    if not fast5_files:
        print(f"  WARNING: No FAST5 files in {fast5_dir}")
        return None

    for f in fast5_files:
        dst = os.path.join(subset_dir, os.path.basename(f))
        if not os.path.exists(dst):
            os.symlink(f, dst)

    print(f"  Prepared subset: {len(fast5_files)} files -> {subset_dir}")
    return subset_dir


def run_rawhash2(dataset_id, segmenter, threads=4):
    """Run RawHash2 indexing + mapping with a specific segmenter."""
    cfg = DATASETS[dataset_id]
    ref_fa = os.path.join(DATA_DIR, dataset_id, "ref.fa")
    fast5_subset = prepare_fast5_subset(dataset_id)

    if fast5_subset is None or not os.path.exists(ref_fa):
        print(f"  SKIP {dataset_id}/{segmenter}: missing data")
        return None

    run_dir = os.path.join(OUTPUT_DIR, dataset_id, segmenter)
    os.makedirs(run_dir, exist_ok=True)

    idx_file = os.path.join(run_dir, "index.ind")
    paf_file = os.path.join(run_dir, "output.paf")

    # Build common args
    seg_args = f"--segmenter {segmenter}" if segmenter != "default" else ""
    extra = cfg["extra"]

    # Index step
    idx_cmd = (f"{RAWHASH2} -x {cfg['preset']} -t {threads} "
               f"-p {cfg['pore']} -d {idx_file} {extra} {seg_args} {ref_fa}")

    print(f"  Indexing {dataset_id} with {segmenter}...")
    t0 = time.time()
    idx_result = subprocess.run(idx_cmd, shell=True, capture_output=True, text=True, timeout=600)
    idx_time = time.time() - t0

    if idx_result.returncode != 0:
        print(f"    Index FAILED: {idx_result.stderr[:200]}")
        # Save error for debugging
        with open(os.path.join(run_dir, "index.err"), 'w') as f:
            f.write(idx_result.stderr)
        return None

    # Map step
    map_cmd = (f"{RAWHASH2} -x {cfg['preset']} -t {threads} "
               f"-o {paf_file} {extra} {seg_args} {idx_file} {fast5_subset}")

    print(f"  Mapping {dataset_id} with {segmenter}...")
    t0 = time.time()
    map_result = subprocess.run(map_cmd, shell=True, capture_output=True, text=True, timeout=1200)
    map_time = time.time() - t0

    if map_result.returncode != 0:
        print(f"    Map FAILED: {map_result.stderr[:200]}")
        with open(os.path.join(run_dir, "map.err"), 'w') as f:
            f.write(map_result.stderr)
        return None

    # Save stderr (contains RawHash2 stats)
    with open(os.path.join(run_dir, "map.log"), 'w') as f:
        f.write(map_result.stderr)

    # Parse PAF output
    metrics = parse_paf(paf_file)
    metrics["index_time"] = idx_time
    metrics["map_time"] = map_time
    metrics["total_time"] = idx_time + map_time

    # Parse RawHash2 stderr for event stats
    stderr_metrics = parse_rawhash2_stderr(map_result.stderr)
    metrics.update(stderr_metrics)

    print(f"    Done: {metrics.get('n_mapped', 0)} mapped, "
          f"{metrics.get('n_unmapped', 0)} unmapped, "
          f"time={metrics['total_time']:.1f}s")

    return metrics


def parse_paf(paf_path):
    """Parse PAF file for mapping statistics."""
    metrics = {
        "n_mapped": 0,
        "n_unmapped": 0,
        "mapq_values": [],
        "alignment_lengths": [],
        "read_lengths": [],
        "mapping_quality_mean": 0,
    }

    if not os.path.exists(paf_path):
        return metrics

    with open(paf_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            read_len = int(parts[1]) if parts[1] != '*' else 0
            metrics["read_lengths"].append(read_len)

            if parts[5] != '*' and parts[11] != '*':
                mapq = int(parts[11])
                aln_len = int(parts[10]) if parts[10] != '*' else 0
                if mapq > 0:
                    metrics["n_mapped"] += 1
                    metrics["mapq_values"].append(mapq)
                    metrics["alignment_lengths"].append(aln_len)
                else:
                    metrics["n_unmapped"] += 1
            else:
                metrics["n_unmapped"] += 1

    if metrics["mapq_values"]:
        metrics["mapping_quality_mean"] = np.mean(metrics["mapq_values"])

    total = metrics["n_mapped"] + metrics["n_unmapped"]
    metrics["mapping_rate"] = metrics["n_mapped"] / total if total > 0 else 0

    return metrics


def parse_rawhash2_stderr(stderr):
    """Extract event detection statistics from RawHash2 stderr."""
    metrics = {}
    for line in stderr.split('\n'):
        line = line.strip()
        if 'mapped' in line.lower() and 'reads' in line.lower():
            # Try to extract mapping counts
            pass
        if 'event' in line.lower():
            # Try to extract event counts
            pass
    return metrics


def load_signals_from_fast5(fast5_dir, max_reads=MAX_READS):
    """Load raw signals from FAST5 files."""
    signals = []
    fast5_files = sorted(glob.glob(os.path.join(fast5_dir, "*.fast5")))[:max_reads]

    for fpath in fast5_files:
        try:
            with h5py.File(fpath, 'r') as f:
                # Try multi-read format
                for key in f.keys():
                    if key.startswith('read_'):
                        sig = f[key]['Raw/Signal'][:].astype(np.float64)
                        signals.append(sig)
                        if len(signals) >= max_reads:
                            return signals
                # Try single-read format
                if 'Raw' in f and 'Reads' in f['Raw']:
                    for read_key in f['Raw/Reads'].keys():
                        sig = f['Raw/Reads'][read_key]['Signal'][:].astype(np.float64)
                        signals.append(sig)
                        if len(signals) >= max_reads:
                            return signals
        except Exception as e:
            continue

    return signals


def compute_segmentation_quality(sig, breakpoints):
    """Compute SNR and R² for a segmentation."""
    if len(breakpoints) == 0 or len(sig) < 10:
        return {"snr": 0, "r2": 0, "n_events": 0, "mean_event_len": len(sig)}

    boundaries = sorted(set([0] + list(breakpoints) + [len(sig)]))
    segments = []
    for i in range(len(boundaries) - 1):
        s, e = boundaries[i], boundaries[i+1]
        if e > s:
            segments.append(sig[s:e])

    if len(segments) < 2:
        return {"snr": 0, "r2": 0, "n_events": len(segments), "mean_event_len": len(sig)}

    # Intra-segment variance (average within-segment variance)
    intra_vars = []
    weighted_intra = 0
    total_n = 0
    for seg in segments:
        if len(seg) > 1:
            v = np.var(seg)
            intra_vars.append(v)
            weighted_intra += v * len(seg)
            total_n += len(seg)

    intra_var = weighted_intra / total_n if total_n > 0 else 1e-10

    # Inter-segment variance (variance of segment means)
    seg_means = [np.mean(seg) for seg in segments]
    inter_var = np.var(seg_means) if len(seg_means) > 1 else 0

    # SNR
    if intra_var > 1e-10:
        snr = 10 * np.log10(inter_var / intra_var)
    else:
        snr = 30.0  # cap

    # R²
    global_var = np.var(sig)
    r2 = 1.0 - (intra_var / global_var) if global_var > 1e-10 else 0

    event_lengths = [len(seg) for seg in segments]

    return {
        "snr": snr,
        "r2": r2,
        "n_events": len(segments),
        "mean_event_len": np.mean(event_lengths),
        "cv_event_len": np.std(event_lengths) / np.mean(event_lengths) if np.mean(event_lengths) > 0 else 0,
    }


def python_segment_signal(sig, method="default"):
    """Segment a signal using Python implementations (for quality metrics)."""
    n = len(sig)
    if n < 20:
        return []

    # MAD normalize
    med = np.median(sig)
    mad = np.median(np.abs(sig - med))
    if mad > 0:
        sig_norm = (sig - med) / (mad * 1.4826)
    else:
        sig_norm = sig - med

    if method == "default":
        return _segment_tstat(sig_norm)
    elif method == "hmm":
        return _segment_hmm(sig_norm)
    elif method == "pelt":
        return _segment_pelt(sig_norm)
    elif method == "binseg":
        return _segment_binseg(sig_norm)
    elif method == "window":
        return _segment_window(sig_norm)
    return []


def _segment_tstat(sig, w1=8, w2=40, t1=1.2, t2=4.5):
    """ONT t-statistic segmenter."""
    n = len(sig)
    ps = np.zeros(n + 1)
    pss = np.zeros(n + 1)
    for i in range(n):
        ps[i+1] = ps[i] + sig[i]
        pss[i+1] = pss[i] + sig[i]**2

    def tstat(ps, pss, n, w):
        eta = np.finfo(np.float32).tiny
        t = np.zeros(n)
        if n < 2*w or w < 2:
            return t
        for i in range(w, n - w + 1):
            s1 = ps[i] - ps[max(0, i-w)]
            sq1 = pss[i] - pss[max(0, i-w)]
            s2 = ps[i+w] - ps[i]
            sq2 = pss[i+w] - pss[i]
            m1, m2 = s1/w, s2/w
            cv = max((sq1/w - m1**2 + sq2/w - m2**2)/w, eta)
            t[i] = abs(m2 - m1) / np.sqrt(cv)
        return t

    t1a = tstat(ps, pss, n, w1)
    breakpoints = []
    for i in range(1, n-1):
        if t1a[i] > t1 and t1a[i] > t1a[i-1] and t1a[i] > t1a[i+1]:
            breakpoints.append(i)
    return breakpoints


def _segment_hmm(sig, n_states=4):
    """HMM Viterbi segmenter."""
    n = len(sig)
    if n < 20:
        return []

    # Initialize with quantiles
    sorted_sig = np.sort(sig)
    state_means = np.array([sorted_sig[(s+1)*n//(n_states+1)] for s in range(n_states)])
    state_stds = np.ones(n_states)

    # K-means (3 iterations)
    assignments = np.zeros(n, dtype=int)
    for _ in range(3):
        for i in range(n):
            dists = np.abs(sig[i] - state_means)
            assignments[i] = np.argmin(dists)
        for s in range(n_states):
            mask = assignments == s
            if mask.sum() > 1:
                state_means[s] = sig[mask].mean()
                state_stds[s] = max(sig[mask].std(), 0.1)

    # Viterbi
    stay_log = np.log(0.95)
    trans_log = np.log(0.05 / (n_states - 1))

    V = np.full((n, n_states), -np.inf)
    path = np.zeros((n, n_states), dtype=int)

    for s in range(n_states):
        emit = -0.5 * np.log(2*np.pi*state_stds[s]**2) - 0.5*((sig[0]-state_means[s])/state_stds[s])**2
        V[0, s] = np.log(1.0/n_states) + emit

    for t in range(1, n):
        for s in range(n_states):
            emit = -0.5 * np.log(2*np.pi*state_stds[s]**2) - 0.5*((sig[t]-state_means[s])/state_stds[s])**2
            for ps in range(n_states):
                v = V[t-1, ps] + (stay_log if ps == s else trans_log)
                if v > V[t, s]:
                    V[t, s] = v + emit
                    path[t, s] = ps

    # Backtrack
    states = np.zeros(n, dtype=int)
    states[-1] = np.argmax(V[-1])
    for t in range(n-2, -1, -1):
        states[t] = path[t+1, states[t+1]]

    return [i for i in range(1, n) if states[i] != states[i-1]]


def _segment_pelt(sig, min_size=5):
    """PELT changepoint detection."""
    n = len(sig)
    if n < 2 * min_size:
        return []

    penalty = 2.0 * np.log(n)
    ps = np.zeros(n + 1)
    pss = np.zeros(n + 1)
    for i in range(n):
        ps[i+1] = ps[i] + sig[i]
        pss[i+1] = pss[i] + sig[i]**2

    F = np.full(n + 1, np.inf)
    F[0] = -penalty
    prev = np.zeros(n + 1, dtype=int)
    cands = [0]

    for j in range(min_size, n + 1):
        new_cands = []
        for t in cands:
            if j - t < min_size:
                continue
            seg_sum = ps[j] - ps[t]
            seg_sumsq = pss[j] - pss[t]
            seg_len = j - t
            seg_mean = seg_sum / seg_len
            cost = seg_sumsq - seg_sum * seg_mean
            total = F[t] + cost + penalty
            if total < F[j]:
                F[j] = total
                prev[j] = t
            if F[t] + cost <= F[j]:
                new_cands.append(t)
        new_cands.append(j)
        cands = new_cands

    # Backtrack
    cps = []
    pos = n
    while pos > 0:
        if prev[pos] > 0:
            cps.append(prev[pos])
        pos = prev[pos]

    return sorted(cps)


def _segment_binseg(sig, min_size=5, max_depth=500):
    """Binary segmentation."""
    n = len(sig)
    if n < 2 * min_size:
        return []

    penalty = 2.0 * np.log(n)
    ps = np.zeros(n + 1)
    pss = np.zeros(n + 1)
    for i in range(n):
        ps[i+1] = ps[i] + sig[i]
        pss[i+1] = pss[i] + sig[i]**2

    cps = []

    def recurse(start, end, depth=0):
        if end - start < 2 * min_size or depth > max_depth or len(cps) > n // min_size:
            return

        total_sum = ps[end] - ps[start]
        total_sumsq = pss[end] - pss[start]
        total_len = end - start
        total_cost = total_sumsq - total_sum**2 / total_len

        best_gain = -1
        best_t = start + min_size

        for t in range(start + min_size, end - min_size + 1):
            left_sum = ps[t] - ps[start]
            left_sumsq = pss[t] - pss[start]
            left_len = t - start
            left_cost = left_sumsq - left_sum**2 / left_len

            right_sum = ps[end] - ps[t]
            right_sumsq = pss[end] - pss[t]
            right_len = end - t
            right_cost = right_sumsq - right_sum**2 / right_len

            gain = total_cost - left_cost - right_cost
            if gain > best_gain:
                best_gain = gain
                best_t = t

        if best_gain > penalty:
            cps.append(best_t)
            recurse(start, best_t, depth + 1)
            recurse(best_t, end, depth + 1)

    recurse(0, n)
    return sorted(cps)


def _segment_window(sig, w=20, min_size=5):
    """Window-based segmenter."""
    n = len(sig)
    if n < 2 * w:
        return []

    penalty = 2.0 * np.log(n)
    ps = np.zeros(n + 1)
    pss = np.zeros(n + 1)
    for i in range(n):
        ps[i+1] = ps[i] + sig[i]
        pss[i+1] = pss[i] + sig[i]**2

    breakpoints = []
    for i in range(w, n - w, min_size):
        left_s = max(0, i - w)
        right_e = min(n, i + w)

        total_sum = ps[right_e] - ps[left_s]
        total_sumsq = pss[right_e] - pss[left_s]
        total_len = right_e - left_s
        total_cost = total_sumsq - total_sum**2 / total_len

        left_sum = ps[i] - ps[left_s]
        left_sumsq = pss[i] - pss[left_s]
        left_len = i - left_s
        left_cost = left_sumsq - left_sum**2 / left_len

        right_sum = ps[right_e] - ps[i]
        right_sumsq = pss[right_e] - pss[i]
        right_len = right_e - i
        right_cost = right_sumsq - right_sum**2 / right_len

        gain = total_cost - left_cost - right_cost
        if gain > penalty:
            breakpoints.append(i)

    return breakpoints


# =============================================================================
# PDF Report Generation
# =============================================================================

def generate_pdf(all_results, seg_quality, pdf_path):
    """Generate comprehensive benchmark PDF."""
    print(f"\nGenerating PDF report: {pdf_path}")

    with PdfPages(pdf_path) as pdf:
        # ---- Page 1: Title ----
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis('off')
        ax.text(0.5, 0.75, "RawHash2 Segmenter Benchmark", fontsize=28,
                fontweight='bold', ha='center', transform=ax.transAxes, family='serif')
        ax.text(0.5, 0.65, "Comparing Event Detection Algorithms for Real-Time Nanopore Mapping",
                fontsize=14, ha='center', transform=ax.transAxes, family='serif', color='#555')

        info_lines = [
            f"Segmenters tested: {', '.join(SEGMENTER_LABELS.values())}",
            f"Datasets: {len(DATASETS)} ({', '.join(d['name'] for d in DATASETS.values())})",
            f"Reads per dataset: {MAX_READS}",
            f"RawHash2 binary: {RAWHASH2}",
            f"Generated: {time.strftime('%Y-%m-%d %H:%M')}",
        ]
        ax.text(0.5, 0.45, '\n'.join(info_lines), fontsize=10, ha='center',
                va='top', transform=ax.transAxes, family='monospace',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f0f0', alpha=0.8))

        # Segmenter categories
        cat_text = "Segmenter Categories:\n"
        for cat, segs in SEGMENTER_CATEGORIES.items():
            cat_text += f"  {cat}: {', '.join(SEGMENTER_LABELS[s] for s in segs)}\n"
        ax.text(0.5, 0.18, cat_text, fontsize=10, ha='center', va='top',
                transform=ax.transAxes, family='monospace')

        ax.text(0.5, 0.03, "P&S Arch4Health - ETH Zurich / SAFARI Research Group",
                fontsize=9, ha='center', transform=ax.transAxes, family='serif', color='#888')
        pdf.savefig(fig); plt.close()

        # ---- Page 2: Mapping Performance Overview ----
        fig = plt.figure(figsize=(16, 10))
        fig.suptitle("Mapping Performance by Segmenter", fontsize=16, fontweight='bold')
        gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

        datasets_with_results = [d for d in DATASETS.keys() if d in all_results and all_results[d]]

        if datasets_with_results:
            # 2a: Mapping rate per segmenter per dataset
            ax = fig.add_subplot(gs[0, 0])
            x = np.arange(len(datasets_with_results))
            width = 0.15
            for i, seg in enumerate(SEGMENTERS):
                rates = []
                for did in datasets_with_results:
                    if seg in all_results.get(did, {}):
                        rates.append(all_results[did][seg].get("mapping_rate", 0) * 100)
                    else:
                        rates.append(0)
                ax.bar(x + i * width, rates, width, label=SEGMENTER_LABELS[seg],
                       color=SEGMENTER_COLORS[seg], alpha=0.8)
            ax.set_ylabel("Mapping Rate (%)", fontsize=10)
            ax.set_title("Mapping Rate", fontsize=12)
            ax.set_xticks(x + width * 2)
            short_names = [DATASETS[d]["name"].split(":")[1].strip()[:12] for d in datasets_with_results]
            ax.set_xticklabels(short_names, fontsize=8, rotation=20)
            ax.legend(fontsize=7, loc='lower left')
            ax.set_ylim(0, 105)
            ax.grid(axis='y', alpha=0.3)

            # 2b: Mean mapping quality
            ax = fig.add_subplot(gs[0, 1])
            for i, seg in enumerate(SEGMENTERS):
                quals = []
                for did in datasets_with_results:
                    if seg in all_results.get(did, {}):
                        quals.append(all_results[did][seg].get("mapping_quality_mean", 0))
                    else:
                        quals.append(0)
                ax.bar(x + i * width, quals, width, label=SEGMENTER_LABELS[seg],
                       color=SEGMENTER_COLORS[seg], alpha=0.8)
            ax.set_ylabel("Mean MAPQ", fontsize=10)
            ax.set_title("Mapping Quality", fontsize=12)
            ax.set_xticks(x + width * 2)
            ax.set_xticklabels(short_names, fontsize=8, rotation=20)
            ax.legend(fontsize=7, loc='lower left')
            ax.grid(axis='y', alpha=0.3)

            # 2c: Execution time
            ax = fig.add_subplot(gs[1, 0])
            for i, seg in enumerate(SEGMENTERS):
                times = []
                for did in datasets_with_results:
                    if seg in all_results.get(did, {}):
                        times.append(all_results[did][seg].get("map_time", 0))
                    else:
                        times.append(0)
                ax.bar(x + i * width, times, width, label=SEGMENTER_LABELS[seg],
                       color=SEGMENTER_COLORS[seg], alpha=0.8)
            ax.set_ylabel("Time (seconds)", fontsize=10)
            ax.set_title("Mapping Time (50 reads)", fontsize=12)
            ax.set_xticks(x + width * 2)
            ax.set_xticklabels(short_names, fontsize=8, rotation=20)
            ax.legend(fontsize=7, loc='upper left')
            ax.grid(axis='y', alpha=0.3)

            # 2d: Number of mapped reads
            ax = fig.add_subplot(gs[1, 1])
            for i, seg in enumerate(SEGMENTERS):
                mapped = []
                for did in datasets_with_results:
                    if seg in all_results.get(did, {}):
                        mapped.append(all_results[did][seg].get("n_mapped", 0))
                    else:
                        mapped.append(0)
                ax.bar(x + i * width, mapped, width, label=SEGMENTER_LABELS[seg],
                       color=SEGMENTER_COLORS[seg], alpha=0.8)
            ax.set_ylabel("Reads Mapped", fontsize=10)
            ax.set_title("Number of Mapped Reads", fontsize=12)
            ax.set_xticks(x + width * 2)
            ax.set_xticklabels(short_names, fontsize=8, rotation=20)
            ax.legend(fontsize=7, loc='lower left')
            ax.grid(axis='y', alpha=0.3)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig); plt.close()

        # ---- Page 3: Segmentation Quality Metrics ----
        if seg_quality:
            fig = plt.figure(figsize=(16, 10))
            fig.suptitle("Segmentation Quality Metrics", fontsize=16, fontweight='bold')
            gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

            datasets_with_sq = [d for d in DATASETS.keys() if d in seg_quality and seg_quality[d]]

            if datasets_with_sq:
                x = np.arange(len(datasets_with_sq))
                short_names_sq = [DATASETS[d]["name"].split(":")[1].strip()[:12] for d in datasets_with_sq]

                # 3a: SNR
                ax = fig.add_subplot(gs[0, 0])
                for i, seg in enumerate(SEGMENTERS):
                    snrs = []
                    for did in datasets_with_sq:
                        if seg in seg_quality.get(did, {}):
                            snrs.append(seg_quality[did][seg].get("snr_mean", 0))
                        else:
                            snrs.append(0)
                    ax.bar(x + i * width, snrs, width, label=SEGMENTER_LABELS[seg],
                           color=SEGMENTER_COLORS[seg], alpha=0.8)
                ax.set_ylabel("SNR (dB)", fontsize=10)
                ax.set_title("Signal-to-Noise Ratio", fontsize=12)
                ax.set_xticks(x + width * 2)
                ax.set_xticklabels(short_names_sq, fontsize=8, rotation=20)
                ax.legend(fontsize=7)
                ax.grid(axis='y', alpha=0.3)

                # 3b: R²
                ax = fig.add_subplot(gs[0, 1])
                for i, seg in enumerate(SEGMENTERS):
                    r2s = []
                    for did in datasets_with_sq:
                        if seg in seg_quality.get(did, {}):
                            r2s.append(seg_quality[did][seg].get("r2_mean", 0))
                        else:
                            r2s.append(0)
                    ax.bar(x + i * width, r2s, width, label=SEGMENTER_LABELS[seg],
                           color=SEGMENTER_COLORS[seg], alpha=0.8)
                ax.set_ylabel("R-squared", fontsize=10)
                ax.set_title("Segmentation R-squared", fontsize=12)
                ax.set_xticks(x + width * 2)
                ax.set_xticklabels(short_names_sq, fontsize=8, rotation=20)
                ax.legend(fontsize=7)
                ax.set_ylim(0, 1.05)
                ax.grid(axis='y', alpha=0.3)

                # 3c: Mean event length
                ax = fig.add_subplot(gs[1, 0])
                for i, seg in enumerate(SEGMENTERS):
                    lens = []
                    for did in datasets_with_sq:
                        if seg in seg_quality.get(did, {}):
                            lens.append(seg_quality[did][seg].get("mean_event_len_mean", 0))
                        else:
                            lens.append(0)
                    ax.bar(x + i * width, lens, width, label=SEGMENTER_LABELS[seg],
                           color=SEGMENTER_COLORS[seg], alpha=0.8)
                ax.set_ylabel("Mean Event Length (samples)", fontsize=10)
                ax.set_title("Average Event Length", fontsize=12)
                ax.set_xticks(x + width * 2)
                ax.set_xticklabels(short_names_sq, fontsize=8, rotation=20)
                ax.legend(fontsize=7)
                ax.grid(axis='y', alpha=0.3)

                # 3d: Number of events
                ax = fig.add_subplot(gs[1, 1])
                for i, seg in enumerate(SEGMENTERS):
                    evts = []
                    for did in datasets_with_sq:
                        if seg in seg_quality.get(did, {}):
                            evts.append(seg_quality[did][seg].get("n_events_mean", 0))
                        else:
                            evts.append(0)
                    ax.bar(x + i * width, evts, width, label=SEGMENTER_LABELS[seg],
                           color=SEGMENTER_COLORS[seg], alpha=0.8)
                ax.set_ylabel("Mean Events per Read", fontsize=10)
                ax.set_title("Number of Events Detected", fontsize=12)
                ax.set_xticks(x + width * 2)
                ax.set_xticklabels(short_names_sq, fontsize=8, rotation=20)
                ax.legend(fontsize=7)
                ax.grid(axis='y', alpha=0.3)

            plt.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig(fig); plt.close()

        # ---- Page 4: Segmentation Quality Boxplots ----
        if seg_quality:
            for metric_key, metric_name, metric_unit in [
                ("snr_all", "SNR", "dB"),
                ("r2_all", "R-squared", ""),
                ("n_events_all", "Events per Read", ""),
            ]:
                fig, axes = plt.subplots(2, 3, figsize=(16, 10))
                fig.suptitle(f"{metric_name} Distribution by Segmenter",
                             fontsize=16, fontweight='bold')

                for idx, did in enumerate(DATASETS.keys()):
                    row, col = idx // 3, idx % 3
                    ax = axes[row][col]

                    if did in seg_quality and seg_quality[did]:
                        data = []
                        labels = []
                        colors = []
                        for seg in SEGMENTERS:
                            if seg in seg_quality[did] and metric_key in seg_quality[did][seg]:
                                vals = seg_quality[did][seg][metric_key]
                                if len(vals) > 0:
                                    data.append(vals)
                                    labels.append(seg[:6])
                                    colors.append(SEGMENTER_COLORS[seg])

                        if data:
                            bp = ax.boxplot(data, patch_artist=True, tick_labels=labels)
                            for patch, color in zip(bp['boxes'], colors):
                                patch.set_facecolor(color)
                                patch.set_alpha(0.6)

                    ax.set_title(DATASETS[did]["name"], fontsize=10)
                    ax.tick_params(labelsize=8)
                    if col == 0:
                        ax.set_ylabel(f"{metric_name} ({metric_unit})" if metric_unit else metric_name, fontsize=9)
                    ax.grid(axis='y', alpha=0.3)

                plt.tight_layout(rect=[0, 0, 1, 0.95])
                pdf.savefig(fig); plt.close()

        # ---- Page 7: Per-dataset detailed results table ----
        for did in DATASETS.keys():
            if did not in all_results or not all_results[did]:
                continue

            fig, ax = plt.subplots(figsize=(14, 8))
            ax.axis('off')
            ax.text(0.5, 0.95, f"Detailed Results: {DATASETS[did]['name']}", fontsize=16,
                    fontweight='bold', ha='center', transform=ax.transAxes)

            # Build table
            headers = ["Segmenter", "Mapped", "Unmapped", "Map Rate", "MAPQ", "Map Time(s)", "Idx Time(s)"]
            rows = []
            for seg in SEGMENTERS:
                if seg in all_results[did]:
                    m = all_results[did][seg]
                    rows.append([
                        SEGMENTER_LABELS[seg],
                        str(m.get("n_mapped", 0)),
                        str(m.get("n_unmapped", 0)),
                        f"{m.get('mapping_rate', 0)*100:.1f}%",
                        f"{m.get('mapping_quality_mean', 0):.1f}",
                        f"{m.get('map_time', 0):.2f}",
                        f"{m.get('index_time', 0):.2f}",
                    ])

            if rows:
                table = ax.table(cellText=rows, colLabels=headers,
                                cellLoc='center', loc='center',
                                colColours=['#4CAF50']*len(headers))
                table.auto_set_font_size(False)
                table.set_fontsize(9)
                table.scale(1.0, 1.5)

                # Color header
                for (i, j), cell in table.get_celld().items():
                    if i == 0:
                        cell.set_text_props(fontweight='bold', color='white')
                        cell.set_facecolor('#2E7D32')
                    elif i % 2 == 0:
                        cell.set_facecolor('#E8F5E9')

            # Add segmentation quality if available
            if did in seg_quality and seg_quality[did]:
                sq_headers = ["Segmenter", "SNR (dB)", "R-squared", "Events", "Evt Len", "CV Evt Len"]
                sq_rows = []
                for seg in SEGMENTERS:
                    if seg in seg_quality[did]:
                        sq = seg_quality[did][seg]
                        sq_rows.append([
                            SEGMENTER_LABELS[seg],
                            f"{sq.get('snr_mean', 0):.2f}",
                            f"{sq.get('r2_mean', 0):.4f}",
                            f"{sq.get('n_events_mean', 0):.0f}",
                            f"{sq.get('mean_event_len_mean', 0):.1f}",
                            f"{sq.get('cv_event_len_mean', 0):.3f}",
                        ])

                if sq_rows:
                    ax2 = fig.add_axes([0.1, 0.05, 0.8, 0.25])
                    ax2.axis('off')
                    ax2.text(0.5, 0.95, "Segmentation Quality", fontsize=12,
                             fontweight='bold', ha='center', transform=ax2.transAxes)
                    table2 = ax2.table(cellText=sq_rows, colLabels=sq_headers,
                                      cellLoc='center', loc='center',
                                      colColours=['#1565C0']*len(sq_headers))
                    table2.auto_set_font_size(False)
                    table2.set_fontsize(9)
                    table2.scale(1.0, 1.5)
                    for (i, j), cell in table2.get_celld().items():
                        if i == 0:
                            cell.set_text_props(fontweight='bold', color='white')
                            cell.set_facecolor('#0D47A1')
                        elif i % 2 == 0:
                            cell.set_facecolor('#E3F2FD')

            pdf.savefig(fig); plt.close()

        # ---- Final page: RawHash2 Benchmarking System Summary ----
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis('off')
        ax.text(0.5, 0.95, "How RawHash2's Benchmarking System Works",
                fontsize=16, fontweight='bold', ha='center', transform=ax.transAxes)

        summary = """
RawHash2 Benchmarking Pipeline:

  1. INDEXING
     rawhash2 -x <preset> -t <threads> -p <pore_model> -d <index> <ref.fa>
     - Builds a hash-table index from the reference genome
     - Index includes quantized k-mer signal patterns from the pore model
     - The --segmenter flag affects how reference signals are segmented into events

  2. MAPPING
     rawhash2 -x <preset> -t <threads> -o <output.paf> <index> <fast5_dir>
     - Reads raw signals from FAST5 files
     - Event detection (segmentation) splits raw signal into discrete events
     - Events are quantized and hashed for seed lookup
     - Seeds are chained via dynamic programming (like minimap2)
     - Output: PAF file with mapping coordinates and quality

  3. GROUND TRUTH (standard benchmark)
     - minimap2 aligns basecalled reads (reads.fasta) to reference -> ground truth PAF
     - RawHash2 PAF is compared to ground truth using compare_pafs.py
     - Metrics: TP (correct map), FP (wrong map), FN (missed), TN (correctly unmapped)
     - Precision = TP/(TP+FP), Recall = TP/(TP+FN), F1 = harmonic mean

  4. TOOLS COMPARED (in standard benchmark)
     - UNCALLED: another signal-level mapper (DTW-based)
     - Sigmap: signal-level mapper (using events + DTW)
     - RawHash / RawHash2: hash-based signal mapper (this tool)
     - minimap2: sequence-level mapper (ground truth)

  This benchmark: We compare SEGMENTERS within RawHash2 itself.
  Same pipeline, same index, different event detection algorithms.
  The segmenter affects both indexing AND mapping quality.
"""
        ax.text(0.05, 0.85, summary, fontsize=9.5, family='monospace',
                va='top', transform=ax.transAxes)
        pdf.savefig(fig); plt.close()

    print(f"PDF saved: {pdf_path}")


# =============================================================================
# Main
# =============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 70)
    print("RawHash2 Segmenter Benchmark")
    print("=" * 70)

    # Check which datasets are available
    available = {}
    for did, cfg in DATASETS.items():
        fast5_dir = os.path.join(DATA_DIR, did, "fast5_files")
        ref_fa = os.path.join(DATA_DIR, did, "ref.fa")
        n_fast5 = len(glob.glob(os.path.join(fast5_dir, "*.fast5")))
        has_ref = os.path.exists(ref_fa)
        available[did] = n_fast5 > 0 and has_ref
        status = f"{n_fast5} fast5, ref={'YES' if has_ref else 'NO'}"
        print(f"  {cfg['name']}: {status} {'[OK]' if available[did] else '[SKIP]'}")

    datasets_to_run = [d for d, ok in available.items() if ok]
    if not datasets_to_run:
        print("\nNo datasets available! Run download_selected.sh first.")
        sys.exit(1)

    print(f"\nRunning benchmarks on {len(datasets_to_run)} datasets with {len(SEGMENTERS)} segmenters")
    print(f"Max reads per dataset: {MAX_READS}")

    # ---- Step 1: Run RawHash2 with each segmenter ----
    all_results = {}
    for did in datasets_to_run:
        print(f"\n{'='*50}")
        print(f"Dataset: {DATASETS[did]['name']}")
        print(f"{'='*50}")

        all_results[did] = {}
        for seg in SEGMENTERS:
            print(f"\n--- Segmenter: {SEGMENTER_LABELS[seg]} ---")
            metrics = run_rawhash2(did, seg, threads=4)
            if metrics:
                all_results[did][seg] = metrics

    # ---- Step 2: Compute segmentation quality ----
    print(f"\n{'='*70}")
    print("Computing Segmentation Quality Metrics")
    print(f"{'='*70}")

    seg_quality = {}
    for did in datasets_to_run:
        print(f"\n  {DATASETS[did]['name']}...")
        fast5_dir = os.path.join(DATA_DIR, did, "fast5_files")
        signals = load_signals_from_fast5(fast5_dir, MAX_READS)

        if not signals:
            print(f"    No signals loaded, skipping quality metrics")
            continue

        print(f"    Loaded {len(signals)} signals")
        seg_quality[did] = {}

        for seg in SEGMENTERS:
            print(f"    {SEGMENTER_LABELS[seg]}...", end=" ", flush=True)
            all_snr, all_r2, all_events, all_evt_len, all_cv = [], [], [], [], []

            t0 = time.time()
            for sig in signals:
                # Subsample very long signals for speed
                if len(sig) > 50000:
                    sig = sig[:50000]

                try:
                    bps = python_segment_signal(sig, seg)
                    q = compute_segmentation_quality(sig, bps)
                    all_snr.append(q["snr"])
                    all_r2.append(q["r2"])
                    all_events.append(q["n_events"])
                    all_evt_len.append(q["mean_event_len"])
                    all_cv.append(q.get("cv_event_len", 0))
                except Exception as e:
                    pass

            elapsed = time.time() - t0

            if all_snr:
                seg_quality[did][seg] = {
                    "snr_mean": np.mean(all_snr),
                    "snr_std": np.std(all_snr),
                    "snr_all": all_snr,
                    "r2_mean": np.mean(all_r2),
                    "r2_std": np.std(all_r2),
                    "r2_all": all_r2,
                    "n_events_mean": np.mean(all_events),
                    "n_events_all": all_events,
                    "mean_event_len_mean": np.mean(all_evt_len),
                    "mean_event_len_all": all_evt_len,
                    "cv_event_len_mean": np.mean(all_cv),
                    "time": elapsed,
                }
                print(f"SNR={np.mean(all_snr):.1f}dB, R2={np.mean(all_r2):.3f}, "
                      f"events={np.mean(all_events):.0f}, time={elapsed:.1f}s")
            else:
                print("no results")

    # ---- Step 3: Save results ----
    results_json = os.path.join(OUTPUT_DIR, "benchmark_results.json")
    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    save_results = {}
    for did in all_results:
        save_results[did] = {}
        for seg in all_results[did]:
            save_results[did][seg] = {k: convert(v) for k, v in all_results[did][seg].items()
                                      if not isinstance(v, list) or len(v) < 100}

    with open(results_json, 'w') as f:
        json.dump(save_results, f, indent=2, default=convert)
    print(f"\nResults saved: {results_json}")

    # ---- Step 4: Generate PDF ----
    pdf_path = os.path.join(OUTPUT_DIR, "segmenter_benchmark.pdf")
    generate_pdf(all_results, seg_quality, pdf_path)

    # Also copy to home
    import shutil
    home_pdf = os.path.expanduser("~/segmenter_benchmark.pdf")
    shutil.copy2(pdf_path, home_pdf)
    print(f"PDF copied to: {home_pdf}")

    print(f"\n{'='*70}")
    print("Benchmark complete!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()

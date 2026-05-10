"""F1 sweep + wall-time benchmark for the Wilcoxon segmenter.

Synthetic benchmark (controlled F1):
  Generate piecewise-constant signal with known boundaries and additive
  Gaussian noise. Sweep (w1, w2, t1, t2, peak_height) and report
  (precision, recall, F1, time_ms) for each cell.

Real-data benchmark (wall-time only):
  Read raw signals from a fast5 file, run the optimal config, report
  ms/read and Msamples/s.
"""
import argparse
import itertools
import json
import os
import sys
import time

import numpy as np

# Allow running from any dir
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from segmenters import (
    tstat_segmenter,
    wilcoxon_segmenter,
)
from boundary_analysis import (
    cross_segmenter_boundary_variance,
    within_segment_sigma2,
)


# -------------------------- synthetic data --------------------------

def gen_synthetic(rng, n=10000, mean_len=9, sigma_level=2.0, sigma_noise=0.5,
                  min_len=3):
    """Piecewise-constant signal + Gaussian noise with known boundaries."""
    parts = []
    bounds = [0]
    while bounds[-1] < n:
        L = max(min_len, int(rng.poisson(mean_len)))
        end = min(bounds[-1] + L, n)
        level = rng.normal(0.0, sigma_level)
        parts.append(rng.normal(level, sigma_noise, end - bounds[-1]))
        bounds.append(end)
    sig = np.concatenate(parts).astype(np.float32)
    return sig, np.asarray(bounds, dtype=np.int64)


# ----------------------------- F1 metric ----------------------------

def f1_boundaries(true_bounds, pred_bounds, tolerance=5):
    """Greedy 1:1 matching of interior boundaries within +-tolerance samples."""
    T = np.asarray(true_bounds[1:-1] if len(true_bounds) > 2 else [], dtype=np.int64)
    P = np.asarray(pred_bounds[1:-1] if len(pred_bounds) > 2 else [], dtype=np.int64)
    if len(T) == 0 and len(P) == 0:
        return 1.0, 0, 0, 0
    if len(P) == 0:
        return 0.0, 0, 0, len(T)
    if len(T) == 0:
        return 0.0, 0, len(P), 0
    matched_p = np.zeros(len(P), dtype=bool)
    tp = 0
    for t in T:
        idx = np.searchsorted(P, t)
        best_j, best_d = -1, tolerance + 1
        for j in (idx - 1, idx):
            if 0 <= j < len(P) and not matched_p[j]:
                d = abs(int(P[j]) - int(t))
                if d <= tolerance and d < best_d:
                    best_j, best_d = j, d
        if best_j >= 0:
            matched_p[best_j] = True
            tp += 1
    fp = int((~matched_p).sum())
    fn = int(len(T) - tp)
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return f1, int(tp), fp, fn


# ----------------------------- sweep -------------------------------

def sweep(sig, true_bounds, grid, tolerance=5):
    results = []
    for w1, w2, t1, t2, ph in itertools.product(*grid):
        if w1 >= w2:
            continue
        t0 = time.perf_counter()
        bounds, _ = wilcoxon_segmenter(sig, w1=w1, w2=w2, t1=t1, t2=t2,
                                       peak_height=ph)
        dt = time.perf_counter() - t0
        f1, tp, fp, fn = f1_boundaries(true_bounds, bounds, tolerance)
        results.append({
            "w1": w1, "w2": w2, "t1": t1, "t2": t2, "peak_height": ph,
            "f1": f1, "tp": tp, "fp": fp, "fn": fn,
            "time_s": dt, "n_pred": len(bounds) - 1,
        })
    return results


# --------------------------- fast5 loader ---------------------------

def load_fast5_signals(path, max_reads=20):
    """Yield (read_id, raw_signal_pA-equivalent-float32) from one fast5."""
    import h5py
    with h5py.File(path, "r") as f:
        keys = list(f.keys())
        for k in keys[:max_reads]:
            grp = f[k]
            if "Raw" in grp:
                rd = grp["Raw"]
            elif "Raw/Reads" in grp:
                rd = list(grp["Raw/Reads"].values())[0]
            else:
                continue
            sig = rd["Signal"][:].astype(np.float32)
            # Center+scale (median MAD) — same as RawHash2 normalization in spirit
            med = np.median(sig)
            mad = np.median(np.abs(sig - med)) * 1.4826 + 1e-6
            sig = (sig - med) / mad
            yield k, sig


# ------------------------------- main -------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results.json")
    ap.add_argument("--n-samples", type=int, default=10000)
    ap.add_argument("--n-trials", type=int, default=5,
                    help="averaged over this many synthetic signals")
    ap.add_argument("--mean-len", type=int, default=9)
    ap.add_argument("--sigma-level", type=float, default=2.0)
    ap.add_argument("--sigma-noise", type=float, default=0.5)
    ap.add_argument("--tolerance", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--fast5", default="",
                    help="optional fast5 for wall-time benchmark on real signal")
    ap.add_argument("--n-real-reads", type=int, default=20)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    grid = (
        [2, 3, 4, 5],                      # w1
        [6, 8, 10, 12, 15],                # w2
        [1.5, 2.0, 2.5, 3.0, 3.5],         # t1
        [1.5, 2.0, 2.5, 3.0],              # t2
        [0.1, 0.2, 0.3, 0.4, 0.5],         # peak_height
    )

    # ------------- 1. Synthetic F1 sweep -------------
    print(f"[sweep] grid size (w1<w2 filtered): "
          f"{sum(1 for w1, w2, *_ in itertools.product(*grid) if w1 < w2)} configs"
          f", trials={args.n_trials}", flush=True)

    per_trial = []
    for trial in range(args.n_trials):
        sig, gt = gen_synthetic(
            rng, n=args.n_samples, mean_len=args.mean_len,
            sigma_level=args.sigma_level, sigma_noise=args.sigma_noise,
        )
        trial_results = sweep(sig, gt, grid, tolerance=args.tolerance)
        per_trial.append(trial_results)
        best = max(trial_results, key=lambda r: r["f1"])
        print(f"[sweep] trial {trial}: best F1 = {best['f1']:.3f} "
              f"@ w1={best['w1']} w2={best['w2']} t1={best['t1']} "
              f"t2={best['t2']} ph={best['peak_height']} "
              f"({best['time_s']*1000:.1f} ms)", flush=True)

    # Aggregate across trials: mean F1 / time per config
    keyfn = lambda r: (r["w1"], r["w2"], r["t1"], r["t2"], r["peak_height"])
    agg = {}
    for trial_results in per_trial:
        for r in trial_results:
            agg.setdefault(keyfn(r), []).append(r)
    sweep_summary = []
    for cfg, runs in agg.items():
        f1s = np.array([r["f1"] for r in runs])
        ts = np.array([r["time_s"] for r in runs])
        sweep_summary.append({
            "w1": cfg[0], "w2": cfg[1], "t1": cfg[2], "t2": cfg[3],
            "peak_height": cfg[4],
            "f1_mean": float(f1s.mean()), "f1_std": float(f1s.std()),
            "time_ms_mean": float(ts.mean() * 1000),
            "n_trials": len(runs),
        })
    sweep_summary.sort(key=lambda r: -r["f1_mean"])
    optimal = sweep_summary[0]
    print(f"[sweep] OPTIMAL (mean over {args.n_trials} trials): "
          f"F1 = {optimal['f1_mean']:.3f} +- {optimal['f1_std']:.3f} "
          f"@ w1={optimal['w1']} w2={optimal['w2']} t1={optimal['t1']} "
          f"t2={optimal['t2']} ph={optimal['peak_height']} "
          f"({optimal['time_ms_mean']:.1f} ms / {args.n_samples} samples)",
          flush=True)

    # ------------- 2. Boundary variance + sigma^2 -------------
    print("[diag] running boundary-variance + sigma^2 across segmenters...",
          flush=True)
    sig, gt = gen_synthetic(
        rng, n=args.n_samples, mean_len=args.mean_len,
        sigma_level=args.sigma_level, sigma_noise=args.sigma_noise,
    )
    seg_outputs = {
        "tstat_default": tstat_segmenter(sig, w1=3, w2=9, t1=4.0, t2=3.5,
                                         peak_height=0.4),
        "tstat_scrappie": tstat_segmenter(sig, w1=3, w2=6, t1=1.4, t2=9.0,
                                          peak_height=0.4),
        "wilcoxon_default": wilcoxon_segmenter(sig, w1=3, w2=9, t1=2.5,
                                               t2=2.0, peak_height=0.3),
        "wilcoxon_optimal": wilcoxon_segmenter(
            sig, w1=optimal["w1"], w2=optimal["w2"],
            t1=optimal["t1"], t2=optimal["t2"],
            peak_height=optimal["peak_height"],
        ),
    }
    diag = {}
    for name, (bounds, _) in seg_outputs.items():
        sigma2_mean, _ = within_segment_sigma2(sig, bounds)
        diag[name] = {
            "n_segments": len(bounds) - 1,
            "noise_sigma2_mean": sigma2_mean,
        }
    pred_list = [seg_outputs[k][0] for k in seg_outputs.keys()]
    _, per_bound, per_seg = cross_segmenter_boundary_variance(
        gt, pred_list, tolerance=20,
    )
    diag["cross_segmenter"] = {
        "names": list(seg_outputs.keys()),
        "per_segmenter": dict(zip(seg_outputs.keys(), per_seg)),
        "per_boundary_std_mean": float(
            np.nanmean(per_bound["std_offset"])
            if per_bound["std_offset"] else float("nan")
        ),
        "per_boundary_std_median": float(
            np.nanmedian(per_bound["std_offset"])
            if per_bound["std_offset"] else float("nan")
        ),
    }
    print("[diag] noise sigma^2 per segmenter:")
    for name in seg_outputs:
        print(f"  {name:20s} sigma^2 = {diag[name]['noise_sigma2_mean']:.4f} "
              f"({diag[name]['n_segments']} segments)")
    print(f"[diag] mean cross-segmenter boundary std = "
          f"{diag['cross_segmenter']['per_boundary_std_mean']:.2f} samples")

    # ------------- 3. Optional real-signal time benchmark -------------
    real_bench = None
    if args.fast5 and os.path.exists(args.fast5):
        print(f"[bench] timing on real reads from {args.fast5}", flush=True)
        per_read_ms = []
        per_read_samples = []
        for rid, sig in load_fast5_signals(args.fast5,
                                           max_reads=args.n_real_reads):
            t0 = time.perf_counter()
            wilcoxon_segmenter(
                sig, w1=optimal["w1"], w2=optimal["w2"],
                t1=optimal["t1"], t2=optimal["t2"],
                peak_height=optimal["peak_height"],
            )
            per_read_ms.append((time.perf_counter() - t0) * 1000)
            per_read_samples.append(len(sig))
        if per_read_ms:
            arr_ms = np.asarray(per_read_ms)
            arr_n = np.asarray(per_read_samples)
            real_bench = {
                "n_reads": len(per_read_ms),
                "ms_per_read_mean": float(arr_ms.mean()),
                "ms_per_read_median": float(np.median(arr_ms)),
                "samples_per_sec": float(arr_n.sum() / (arr_ms.sum() / 1000.0)),
                "Msamples_per_sec": float(
                    arr_n.sum() / (arr_ms.sum() / 1000.0) / 1e6
                ),
            }
            print(f"[bench] {real_bench['n_reads']} reads, "
                  f"{real_bench['ms_per_read_mean']:.1f} ms/read mean "
                  f"({real_bench['Msamples_per_sec']:.2f} Msamples/s)")

    # ------------- 4. Dump JSON -------------
    payload = {
        "config": vars(args),
        "optimal_wilcoxon": optimal,
        "sweep_top20": sweep_summary[:20],
        "diagnostics": diag,
        "real_signal_benchmark": real_bench,
    }
    with open(args.out, "w") as f:
        json.dump(payload, f, indent=2, default=float)
    print(f"[done] results -> {args.out}")


if __name__ == "__main__":
    main()

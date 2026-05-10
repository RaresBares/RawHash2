"""Bench Wilcoxon (default-aligned) and t-stat on 10kb and 100kb signals.

Two regimes per length:
  (a) SYNTHETIC step+Gaussian noise (known boundaries) -> F1 + time
  (b) REAL nanopore signal from fast5 (no GT) -> time + n_segments + sigma^2

Lengths tested:
  10kb  bases  ~  90,000 samples (R9.4 ~9 samples/base)
  100kb bases ~ 900,000 samples
"""
import sys, os, time, glob
import numpy as np
import h5py

sys.path.insert(0, "/tmp/wilcox_seg")
from segmenters import tstat_segmenter, wilcoxon_segmenter
from sweep_f1 import gen_synthetic, f1_boundaries
from boundary_analysis import within_segment_sigma2


def time_segmenter(name, fn, sig, n_runs=3, **kw):
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        bounds, _ = fn(sig, **kw)
        times.append(time.perf_counter() - t0)
    return min(times), bounds


def bench_synth(n_samples, mean_len=9, sigma_noise=0.5, n_trials=5):
    rng = np.random.default_rng(42)
    rows = []
    for trial in range(n_trials):
        sig, gt = gen_synthetic(rng, n=n_samples, mean_len=mean_len,
                                sigma_level=2.0, sigma_noise=sigma_noise)
        # default t-stat
        t_t, bt = time_segmenter("tstat", tstat_segmenter, sig,
                                 w1=3, w2=9, t1=4.0, t2=3.5, peak_height=0.4)
        f1_t, *_ = f1_boundaries(gt, bt, tolerance=5)
        # default-aligned wilcoxon
        t_w, bw = time_segmenter("wilcox", wilcoxon_segmenter, sig,
                                 w1=3, w2=9, t1=1.8, t2=3.0, peak_height=0.4)
        f1_w, *_ = f1_boundaries(gt, bw, tolerance=5)
        rows.append({
            "n_samples": n_samples,
            "trial": trial,
            "tstat_time_s": t_t,
            "tstat_n_seg": len(bt) - 1,
            "tstat_f1": f1_t,
            "wilcox_time_s": t_w,
            "wilcox_n_seg": len(bw) - 1,
            "wilcox_f1": f1_w,
        })
    return rows


def get_real_reads(fast5_paths, target_n, n_reads=5, tolerance=0.3):
    """Pick reads with ~target_n samples (within +-tolerance fraction)."""
    candidates = []
    lo, hi = target_n * (1 - tolerance), target_n * (1 + tolerance)
    for fp in fast5_paths:
        try:
            with h5py.File(fp, "r") as f:
                if "Raw/Reads" in f:
                    for rk in f["Raw/Reads"]:
                        sig = f[f"Raw/Reads/{rk}/Signal"][:]
                        if lo <= len(sig) <= hi:
                            candidates.append((len(sig), sig.astype(np.float32)))
                            if len(candidates) >= n_reads:
                                return candidates
        except Exception:
            pass
    return candidates


def normalize_mad(sig):
    med = np.median(sig)
    mad = np.median(np.abs(sig - med)) * 1.4826 + 1e-9
    return (sig - med) / mad


def bench_real(reads, label):
    rows = []
    for i, (n_orig, sig) in enumerate(reads):
        sig_n = normalize_mad(sig)
        t_t, bt = time_segmenter("tstat", tstat_segmenter, sig_n,
                                 w1=3, w2=9, t1=4.0, t2=3.5, peak_height=0.4)
        s2_t, _ = within_segment_sigma2(sig, bt)  # raw pA-like scale
        t_w, bw = time_segmenter("wilcox", wilcoxon_segmenter, sig_n,
                                 w1=3, w2=9, t1=1.8, t2=3.0, peak_height=0.4)
        s2_w, _ = within_segment_sigma2(sig, bw)
        rows.append({
            "label": label,
            "read": i,
            "n_samples": n_orig,
            "tstat_time_s": t_t, "tstat_n_seg": len(bt) - 1,
            "tstat_sigma2": float(s2_t),
            "wilcox_time_s": t_w, "wilcox_n_seg": len(bw) - 1,
            "wilcox_sigma2": float(s2_w),
        })
    return rows


def fmt_rows(label, rows, synth=True):
    if not rows: return
    arr = lambda k: np.array([r[k] for r in rows])
    print(f"=== {label}  (n={len(rows)}) ===")
    if synth:
        print(f"  tstat:    F1 = {arr('tstat_f1').mean():.3f} +- {arr('tstat_f1').std():.3f}   "
              f"time = {arr('tstat_time_s').mean()*1000:7.2f} ms   "
              f"n_seg = {arr('tstat_n_seg').mean():.0f}")
        print(f"  wilcox:   F1 = {arr('wilcox_f1').mean():.3f} +- {arr('wilcox_f1').std():.3f}   "
              f"time = {arr('wilcox_time_s').mean()*1000:7.2f} ms   "
              f"n_seg = {arr('wilcox_n_seg').mean():.0f}")
        print(f"  slowdown wilcox/tstat = {arr('wilcox_time_s').mean()/arr('tstat_time_s').mean():.1f}x")
    else:
        print(f"  median read len: {int(np.median(arr('n_samples')))} samples")
        print(f"  tstat:    time = {arr('tstat_time_s').mean()*1000:7.2f} ms   "
              f"n_seg = {arr('tstat_n_seg').mean():.0f}   "
              f"sigma2 = {arr('tstat_sigma2').mean():.2f}")
        print(f"  wilcox:   time = {arr('wilcox_time_s').mean()*1000:7.2f} ms   "
              f"n_seg = {arr('wilcox_n_seg').mean():.0f}   "
              f"sigma2 = {arr('wilcox_sigma2').mean():.2f}")
        print(f"  slowdown wilcox/tstat = {arr('wilcox_time_s').mean()/arr('tstat_time_s').mean():.1f}x")
    print()


if __name__ == "__main__":
    print("Bench: 10kb (~90k samples) and 100kb (~900k samples) signals")
    print()

    # ---- SYNTH ----
    rows10k = bench_synth(90000, n_trials=5)
    rows100k = bench_synth(900000, n_trials=3)  # fewer trials, more expensive
    fmt_rows("SYNTH 10kb (90k samples)", rows10k, synth=True)
    fmt_rows("SYNTH 100kb (900k samples)", rows100k, synth=True)

    # ---- REAL ----
    print("Loading real reads...")
    d3_files = sorted(glob.glob("/tmp/d3_fast5/*.fast5"))
    d4_files = sorted(glob.glob("/tmp/d4_fast5/*.fast5"))
    reads_10k = get_real_reads(d3_files + d4_files, 90000, n_reads=5)
    reads_100k = get_real_reads(d4_files, 750000, n_reads=3, tolerance=0.5)
    print(f"  10kb-target reads found: {len(reads_10k)}")
    print(f"  100kb-target reads found: {len(reads_100k)} "
          f"(actual lens: {[r[0] for r in reads_100k]})")
    print()

    if reads_10k:
        fmt_rows("REAL 10kb (~90k samples)", bench_real(reads_10k, "10kb"), synth=False)
    if reads_100k:
        fmt_rows("REAL 100kb (~750k actual)", bench_real(reads_100k, "100kb"), synth=False)

"""Boundary-placement variance + within-segment noise sigma^2.

Two diagnostics:

1. cross_segmenter_boundary_variance: for each ground-truth boundary,
   find the nearest predicted boundary from each segmenter (within a
   tolerance) and report the std-dev of those offsets across segmenters.
   Tells you where segmenters DISAGREE on the boundary location.

2. within_segment_sigma2: for a given (signal, boundary list), the mean
   per-segment residual variance (signal - segment_mean)^2. This is the
   "noise level" each segmenter implicitly assumes when it draws cuts.
"""
import numpy as np


def within_segment_sigma2(sig, bounds, min_seg_len=2):
    """Mean within-segment variance and per-segment array."""
    sig = np.asarray(sig, dtype=np.float64)
    sigsq = []
    for i in range(len(bounds) - 1):
        a, b = int(bounds[i]), int(bounds[i + 1])
        if b - a >= min_seg_len:
            sigsq.append(sig[a:b].var(ddof=1))
    if not sigsq:
        return float("nan"), np.array([], dtype=np.float64)
    arr = np.asarray(sigsq, dtype=np.float64)
    return float(arr.mean()), arr


def cross_segmenter_boundary_variance(reference_bounds, predicted_bounds_list,
                                      tolerance=20):
    """For each interior reference boundary, gather the nearest predicted
    boundary from every segmenter (offset signed). Returns:

        offsets    : (n_ref, n_seg) array, NaN where no match within tol
        per_bound  : dict with std_offset, mean_abs_offset per ref boundary
        per_seg    : dict with mean_offset, std_offset, recall per segmenter
    """
    R = np.asarray(reference_bounds, dtype=np.int64)
    if len(R) > 2:
        R_int = R[1:-1]
    else:
        R_int = R
    n_ref = len(R_int)
    n_seg = len(predicted_bounds_list)
    offsets = np.full((n_ref, n_seg), np.nan, dtype=np.float64)
    for k, P in enumerate(predicted_bounds_list):
        if len(P) == 0:
            continue
        P = np.asarray(P, dtype=np.int64)
        if len(P) > 2:
            P_int = P[1:-1]
        else:
            P_int = P
        if len(P_int) == 0:
            continue
        for j, r in enumerate(R_int):
            idx = np.searchsorted(P_int, r)
            best = None
            for cand in (idx - 1, idx):
                if 0 <= cand < len(P_int):
                    d = int(P_int[cand]) - int(r)
                    if best is None or abs(d) < abs(best):
                        best = d
            if best is not None and abs(best) <= tolerance:
                offsets[j, k] = best
    per_bound = {
        "std_offset": np.nanstd(offsets, axis=1, ddof=0).tolist(),
        "mean_abs_offset": np.nanmean(np.abs(offsets), axis=1).tolist(),
        "n_segmenters_matched": (~np.isnan(offsets)).sum(axis=1).tolist(),
    }
    per_seg = []
    for k in range(n_seg):
        col = offsets[:, k]
        valid = ~np.isnan(col)
        per_seg.append({
            "mean_offset": float(np.nanmean(col)) if valid.any() else float("nan"),
            "std_offset": float(np.nanstd(col, ddof=0)) if valid.any() else float("nan"),
            "recall_within_tol": float(valid.mean()),
        })
    return offsets, per_bound, per_seg

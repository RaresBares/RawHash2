"""Sliding-window event segmenters for nanopore raw signal.

Two implementations sharing the dual-window peak-detection scaffold from
RawHash2's revent.c:
  - tstat_segmenter   : matches the C default (Welch-style t-statistic)
  - wilcoxon_segmenter: replaces t-stat with a Mann-Whitney U |z|-score

Both feed |statistic| into a threshold-crossing peak finder over two
window sizes (short + long), then turn peaks into segment means.

The Wilcoxon variant uses scipy.stats.mannwhitneyu over a sliding view
(asymptotic mode + tie correction), so it is exact w.r.t. the standard
asymptotic Mann-Whitney test.
"""
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.stats import mannwhitneyu, norm


# ------------------------------- statistics ------------------------------

def _prefix_sums(sig):
    s = sig.astype(np.float64, copy=False)
    ps = np.empty(len(s) + 1, dtype=np.float64)
    pss = np.empty(len(s) + 1, dtype=np.float64)
    ps[0] = 0.0
    pss[0] = 0.0
    np.cumsum(s, out=ps[1:])
    np.cumsum(s * s, out=pss[1:])
    return ps, pss


def comp_tstat(sig, w):
    """Welch t-statistic from prefix sums (mirrors revent.c::comp_tstat).

    Returns array of length n+1 with |t| at index i = boundary between
    [i-w, i] and [i, i+w]. Zeros at borders.
    """
    n = len(sig)
    out = np.zeros(n + 1, dtype=np.float64)
    if n < 2 * w or w < 2:
        return out
    ps, pss = _prefix_sums(sig)
    i = np.arange(w, n - w + 1)
    sum1 = ps[i] - ps[i - w]
    ssq1 = pss[i] - pss[i - w]
    sum2 = ps[i + w] - ps[i]
    ssq2 = pss[i + w] - pss[i]
    mean1 = sum1 / w
    mean2 = sum2 / w
    var1 = np.maximum(ssq1 / w - mean1 * mean1, 1e-10)
    var2 = np.maximum(ssq2 / w - mean2 * mean2, 1e-10)
    combined = (var1 + var2) / w
    out[i] = np.abs(mean2 - mean1) / np.sqrt(combined)
    return out


def comp_wilcoxon_z(sig, w):
    """|z|-score of Mann-Whitney U comparing [i-w, i] vs [i, i+w].

    Uses scipy.stats.mannwhitneyu with axis=1 over a sliding view; |z|
    is recovered from the two-sided p-value. Cost dominated by C-side
    rank computation in scipy: ~O(n * w log w).
    """
    n = len(sig)
    out = np.zeros(n + 1, dtype=np.float64)
    if n < 2 * w or w < 2:
        return out
    s = sig.astype(np.float64, copy=False)
    win = sliding_window_view(s, 2 * w)
    left = win[:, :w]
    right = win[:, w:]
    res = mannwhitneyu(left, right, axis=1, alternative="two-sided",
                       method="asymptotic")
    p = np.clip(res.pvalue, 1e-300, 1.0)
    z = np.abs(norm.ppf(p / 2.0))
    out[w:n - w + 1] = z
    return out


# --------------------------- peak detection -----------------------------

def gen_peaks(stat_arrs, thresholds, window_lengths, peak_height):
    """Multi-window threshold-crossing peak finder (port of revent.c::gen_peaks).

    State machine: each detector first descends into a valley, then waits
    for an above-threshold rise of >= peak_height before emitting the peak.
    A confirmed peak masks the other detectors for window_length samples.
    """
    n = len(stat_arrs[0])
    n_d = len(stat_arrs)
    DEF_POS = -1
    DEF_VAL = float("inf")
    state = [
        {
            "sig": stat_arrs[k],
            "thr": float(thresholds[k]),
            "wlen": int(window_lengths[k]),
            "masked_to": 0,
            "peak_pos": DEF_POS,
            "peak_val": DEF_VAL,
            "valid": False,
        }
        for k in range(n_d)
    ]
    peaks = []
    short_wlen = state[0]["wlen"]
    for i in range(n):
        for k in range(n_d):
            d = state[k]
            if i < d["masked_to"]:
                continue
            cur = d["sig"][i]
            if d["peak_pos"] == DEF_POS:
                if cur < d["peak_val"]:
                    d["peak_val"] = cur
                elif cur - d["peak_val"] > peak_height:
                    d["peak_val"] = cur
                    d["peak_pos"] = i
            else:
                if cur > d["peak_val"]:
                    d["peak_val"] = cur
                    d["peak_pos"] = i
                if d["peak_val"] > d["thr"]:
                    for k2 in range(n_d):
                        if k2 != k:
                            state[k2]["masked_to"] = d["peak_pos"] + short_wlen
                            state[k2]["peak_pos"] = DEF_POS
                            state[k2]["peak_val"] = DEF_VAL
                            state[k2]["valid"] = False
                if (
                    d["peak_val"] - cur > peak_height
                    and d["peak_val"] > d["thr"]
                ):
                    d["valid"] = True
                if d["valid"] and (i - d["peak_pos"]) > d["wlen"] // 2:
                    peaks.append(d["peak_pos"])
                    d["peak_pos"] = DEF_POS
                    d["peak_val"] = cur
                    d["valid"] = False
    peaks.sort()
    return np.asarray(peaks, dtype=np.int64)


def gen_events(sig, peaks):
    """Per-segment mean (port of revent.c::gen_events). Returns (bounds, events)."""
    n = len(sig)
    if len(peaks) == 0:
        bounds = np.array([0, n], dtype=np.int64)
    else:
        bounds = np.concatenate([[0], peaks, [n]]).astype(np.int64)
    events = np.empty(len(bounds) - 1, dtype=np.float64)
    for i in range(len(events)):
        a, b = bounds[i], bounds[i + 1]
        events[i] = sig[a:b].mean() if b > a else 0.0
    return bounds, events


# --------------------------- public segmenters ---------------------------

def tstat_segmenter(sig, w1=3, w2=9, t1=4.0, t2=3.5, peak_height=0.4):
    """Default RawHash2 segmenter (Python reference port)."""
    s1 = comp_tstat(sig, w1)
    s2 = comp_tstat(sig, w2)
    peaks = gen_peaks([s1, s2], [t1, t2], [w1, w2], peak_height)
    return gen_events(sig, peaks)


def wilcoxon_segmenter(sig, w1=3, w2=9, t1=2.5, t2=2.0, peak_height=0.3):
    """Sliding-window Mann-Whitney U segmenter."""
    s1 = comp_wilcoxon_z(sig, w1)
    s2 = comp_wilcoxon_z(sig, w2)
    peaks = gen_peaks([s1, s2], [t1, t2], [w1, w2], peak_height)
    return gen_events(sig, peaks)

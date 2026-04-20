#!/usr/bin/env python3
"""Ground-truth evaluation using RawHash2's official method (pafstats.py).

A pair = (read_id, reference_id). Classification:
- TP: pair in both query PAF and truth PAF
- FP: pair only in query (query claims mapping that truth doesn't)
- FN: pair only in truth (query missed what truth found)
- TN: pair in neither

Usage: ground_truth_rawhash_style.py <query.paf> <truth.paf> [subset_reads.txt]
"""
import sys, os, json

def read_paf_pairs(path):
    """Return (mapped_pairs_set, unmapped_pairs_set)."""
    mapped = set(); unmapped = set()
    if not os.path.exists(path):
        return mapped, unmapped
    with open(path) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 6:
                continue
            q = p[0]; r = p[5]
            if r != '*':
                mapped.add((q, r))
            else:
                unmapped.add((q, r))
    return mapped, unmapped


def classify(query_paf, truth_paf, subset_reads=None):
    q_mapped, q_unmapped = read_paf_pairs(query_paf)
    t_mapped, t_unmapped = read_paf_pairs(truth_paf)

    # Optionally restrict to a subset of reads
    if subset_reads is not None:
        q_mapped = {p for p in q_mapped if p[0] in subset_reads}
        q_unmapped = {p for p in q_unmapped if p[0] in subset_reads}
        t_mapped = {p for p in t_mapped if p[0] in subset_reads}
        t_unmapped = {p for p in t_unmapped if p[0] in subset_reads}

    counts = {'tp': 0, 'fp': 0, 'fn': 0, 'tn': 0}
    all_pairs = q_mapped | q_unmapped | t_mapped | t_unmapped
    for pair in all_pairs:
        if pair in q_mapped:
            if pair in t_mapped:
                counts['tp'] += 1
            else:
                counts['fp'] += 1
        elif pair in t_mapped:
            counts['fn'] += 1
        else:
            counts['tn'] += 1

    tp, fp, fn, tn = counts['tp'], counts['fp'], counts['fn'], counts['tn']
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2*tp / (2*tp + fp + fn) if (2*tp + fp + fn) > 0 else 0.0

    return {
        'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn,
        'precision': round(precision, 4),
        'recall': round(recall, 4),
        'f1': round(f1, 4),
        'n_query_mapped': len(q_mapped),
        'n_truth_mapped': len(t_mapped),
        'method': 'rawhash2_pafstats (read_id, ref_id) pair classification',
    }


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: ground_truth_rawhash_style.py <query.paf> <truth.paf> [subset_reads.txt]")
        sys.exit(1)
    q = sys.argv[1]; t = sys.argv[2]
    subset = None
    if len(sys.argv) > 3 and sys.argv[3] and os.path.exists(sys.argv[3]):
        with open(sys.argv[3]) as f:
            subset = set(l.strip() for l in f if l.strip())
    res = classify(q, t, subset)
    print(json.dumps(res, indent=2))

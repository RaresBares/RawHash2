#!/usr/bin/env python3
"""Compare RawHash2 PAF output against minimap2 ground-truth PAF.

Usage: ground_truth_eval.py <truth.paf> <query.paf> [--overlap-frac 0.1]

Returns JSON with TP/FP/FN, Precision, Recall, F1.

Based on the approach in RawHash2's test/scripts/compare_pafs.sh:
- A query mapping is TRUE POSITIVE if it's on the same chromosome/strand as truth
  AND overlaps by at least `overlap_frac` of truth length.
- A truth mapping with no matching query mapping = FALSE NEGATIVE
- A query mapping not matching any truth = FALSE POSITIVE
"""
import sys
import os
import json
from collections import defaultdict

def parse_paf(path, min_mapq=1):
    """Returns dict: read_id -> list of (target, strand, tstart, tend, mapq)"""
    reads = defaultdict(list)
    if not path or path == "-":
        return reads
    with open(path) as f:
        for line in f:
            c = line.strip().split('\t')
            if len(c) < 12:
                continue
            read = c[0]
            if c[5] == '*':
                continue
            try:
                target = c[5]
                strand = c[4]
                tstart = int(c[7])
                tend = int(c[8])
                mapq = int(c[11])
            except (ValueError, IndexError):
                continue
            if mapq < min_mapq:
                continue
            reads[read].append((target, strand, tstart, tend, mapq))
    return reads

def best_mapping(mappings):
    """Return highest-MAPQ mapping (or longest if tie)."""
    return max(mappings, key=lambda m: (m[4], m[3]-m[2]))

def overlap_frac(truth, query):
    """Fraction of TRUTH interval covered by QUERY."""
    ttarget, tstrand, ts, te = truth[:4]
    qtarget, qstrand, qs, qe = query[:4]
    if ttarget != qtarget:
        return 0.0
    if tstrand != qstrand:
        return 0.0
    ovl = max(0, min(te, qe) - max(ts, qs))
    tlen = te - ts
    if tlen == 0:
        return 0.0
    return ovl / tlen


def evaluate(truth_paf, query_paf, overlap_threshold=0.1, subset_reads=None):
    truth = parse_paf(truth_paf, min_mapq=1)
    query = parse_paf(query_paf, min_mapq=1)

    # Take best mapping per read
    truth_best = {r: best_mapping(ms) for r, ms in truth.items()}
    query_best = {r: best_mapping(ms) for r, ms in query.items()}

    # If subset specified, restrict evaluation to reads actually benchmarked
    if subset_reads is not None:
        truth_best = {r: m for r, m in truth_best.items() if r in subset_reads}
        query_best = {r: m for r, m in query_best.items() if r in subset_reads}

    all_reads = set(truth_best) | set(query_best)

    TP, FP, FN = 0, 0, 0
    correct_pos, wrong_pos = 0, 0
    mapq_correct = []
    mapq_wrong = []

    for read in all_reads:
        t = truth_best.get(read)
        q = query_best.get(read)
        if t and q:
            # Both mapped. Check overlap.
            ovl = overlap_frac(t, q)
            if ovl >= overlap_threshold:
                TP += 1
                correct_pos += 1
                mapq_correct.append(q[4])
            else:
                FP += 1  # query wrong position
                FN += 1  # truth missed by query (wrong call)
                wrong_pos += 1
                mapq_wrong.append(q[4])
        elif t and not q:
            FN += 1  # truth exists but query missed
        elif q and not t:
            FP += 1  # query has mapping but truth doesn't

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        'n_truth_reads': len(truth_best),
        'n_query_reads': len(query_best),
        'n_both_mapped': correct_pos + wrong_pos,
        'TP': TP,
        'FP': FP,
        'FN': FN,
        'correct_position': correct_pos,
        'wrong_position': wrong_pos,
        'precision': round(precision, 4),
        'recall': round(recall, 4),
        'f1': round(f1, 4),
        'mean_mapq_correct': round(sum(mapq_correct)/len(mapq_correct), 2) if mapq_correct else 0,
        'mean_mapq_wrong': round(sum(mapq_wrong)/len(mapq_wrong), 2) if mapq_wrong else 0,
        'overlap_threshold': overlap_threshold,
    }


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: ground_truth_eval.py <truth.paf> <query.paf> [overlap_frac] [subset_reads.txt]")
        sys.exit(1)
    truth = sys.argv[1]
    query = sys.argv[2]
    ovl = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
    subset = None
    if len(sys.argv) > 4 and sys.argv[4] and os.path.exists(sys.argv[4]):
        with open(sys.argv[4]) as f:
            subset = set(l.strip() for l in f if l.strip())
    result = evaluate(truth, query, ovl, subset)
    print(json.dumps(result, indent=2))

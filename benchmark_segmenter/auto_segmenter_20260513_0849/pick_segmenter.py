#!/usr/bin/env python3
"""
Rule-based segmenter selector for rawhash2.

Given a dataset's pore chemistry and reference size (and optionally GC%),
returns the (segmenter, mode) pair that the v4 clean-bench data shows is
either the accuracy oracle or Pareto-optimal vs the default segmenter.

Mode is "static" or "auto":
- "auto" sets RH2_TOPK_AUTO=1, enabling the K-auto bootstrap: rawhash2 runs
  the default segmenter on the first 5% of each read to estimate event
  density K, then runs the chosen top-K segmenter with that K.
- "auto" only changes behavior for top-K segmenters (binseg, mad, window,
  gradient, scrappie). For HMM / PELT / BOCD / default it is a no-op.

The selector intentionally leaves RH2_TOPK_AUTO_BOOST at its default (1.0).
The source comment in revent.c recommends BOOST=1.2-1.5 for R10.4 chemistries,
but an overnight sweep on D6 (R10.4), RB_ecoli/RB_dmel/RB_hsap_sub (R10.4.1)
showed that BOOST > 1.0 is *strictly* harmful: accuracy drops by 5-30 pt and
wall time grows 1.2-2.7x. See auto_segmenter_v5_report.pdf, "BOOST effect"
page, for the data.

Rule derived from seq_benchmark/segmenter_paper_bench/v4_bench/v4_aggregated.json
plus the overnight night_runs/ extension (eight datasets, clean per-segmenter
wall times on exclusive nodes, n>=2 reps on D3/D4/RB_ecoli).
"""
from __future__ import annotations
import argparse
import os
import sys
from pathlib import Path


def fasta_stats(ref_path: Path) -> tuple[float, float]:
    """Return (ref_mb, gc_percent) by streaming the FASTA. No biopython dep."""
    n_total = 0
    n_gc = 0
    n_at = 0
    with open(ref_path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                continue
            s = line.rstrip()
            n_total += len(s)
            for b in s:
                if b in b"GgCc":
                    n_gc += 1
                elif b in b"AaTtUu":
                    n_at += 1
    ref_mb = n_total / 1_000_000
    gc_pct = 100.0 * n_gc / (n_gc + n_at) if (n_gc + n_at) else 0.0
    return ref_mb, gc_pct


def pick(chemistry: str, ref_mb: float, gc_pct: float | None = None) -> tuple[str, str]:
    """
    Return (segmenter_name, mode) for the given dataset characteristics.

    chemistry: "R9.4", "R10.4", or "R10.4.1".
    ref_mb:    reference size in megabases.
    gc_pct:    GC content of the reference (0-100). Currently unused but kept
               in the signature for future extension.
    """
    chem = chemistry.strip().upper().replace("R", "R")  # normalize
    if chem == "R9.4":
        # D1 (SARS-CoV-2, 30 kb): default is best by 2 pt and ties on speed.
        # D2/D3/D4 (bacterial, yeast, algae): HMM wins accuracy by 13-26 pt
        # at essentially the same wall time as default.
        if ref_mb < 0.1:
            return ("default", "static")
        return ("hmm", "static")
    if chem == "R10.4":
        # D6 (E. coli R10.4, 5.3 Mb ref, 98 GB fast5): binseg+Auto-K matches
        # PELT's accuracy (within 0.1 pt) while running 9x faster than PELT
        # and marginally faster than default.
        # D7 (S. aureus R10.4, 2.8 Mb ref, low GC 32.8%): default outperforms
        # every alternative segmenter by 3-10 pt - non-default segmenters all
        # land at ovl01 < 0.1. Split at ref_mb = 4 keeps both decisions clean.
        if ref_mb < 4:
            return ("default", "static")
        return ("binseg", "auto")
    if chem == "R10.4.1":
        if ref_mb < 10:
            # RB_ecoli: HMM gives +4.5 pt over default at the same speed.
            return ("hmm", "static")
        if ref_mb < 1000:
            # RB_dmel (143 MB): window with Auto-K beats default by 3 pt
            # and is 14% faster.
            return ("window", "auto")
        # Whole-genome eukaryote (RB_hsap_full, 3.1 Gb): default is robust.
        return ("default", "static")
    return ("default", "static")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--ref", type=Path, required=True, help="Reference FASTA")
    ap.add_argument("--chemistry", required=True,
                    choices=["R9.4", "R10.4", "R10.4.1"],
                    help="Pore chemistry tag")
    ap.add_argument("--shell", action="store_true",
                    help="Emit shell exports instead of human text")
    args = ap.parse_args(argv)

    if not args.ref.exists():
        print(f"ref not found: {args.ref}", file=sys.stderr)
        return 2
    ref_mb, gc = fasta_stats(args.ref)
    seg, mode = pick(args.chemistry, ref_mb, gc)

    if args.shell:
        # Shell-eval'able exports for use inside a sbatch script.
        seg_flag = "" if seg == "default" else f"--segmenter {seg}"
        topk_env = "1" if mode == "auto" else "0"
        print(f"export AUTO_SEG='{seg}'")
        print(f"export AUTO_MODE='{mode}'")
        print(f"export AUTO_SEG_FLAG='{seg_flag}'")
        print(f"export RH2_TOPK_AUTO='{topk_env}'")
        print(f"export AUTO_REF_MB='{ref_mb:.3f}'")
        print(f"export AUTO_REF_GC='{gc:.2f}'")
    else:
        print(f"ref_mb={ref_mb:.3f}  gc={gc:.2f}  chemistry={args.chemistry}")
        print(f"picked: segmenter={seg}  mode={mode}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

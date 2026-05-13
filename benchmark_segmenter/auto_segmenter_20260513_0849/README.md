# Auto-segmenter benchmark 20260513_0849

Rule v2 + Auto-K BOOST sweep + reproducibility run.

## Files

- `auto_segmenter_v5_report.pdf` - paper-grade graphs (no narrative).
- `auto_segmenter_v5_summary.csv` - per-dataset table: default vs rule
  vs oracle, with median wall, accuracy delta, speedup, and run-to-run
  spread.
- `pick_segmenter.py` - the rule (chemistry + ref size -> segmenter, mode).
- `rh2_auto.sh` - wrapper that applies the rule and invokes rawhash2.
- `aggregate_night.py`, `build_report_v2.py` - reproducible pipeline
  used to build the PDF from the raw v4 + overnight CSVs.

## Bench protocol

Every wall-time measurement is from a single sbatch on an exclusive
`cpu_part` node (8 CPUs, 32 GB), exactly the same as
`seq_benchmark/scripts/gt_pipeline/bench_v4_full.sh`.  Multiple repeats
are reported as median + min/max whiskers.

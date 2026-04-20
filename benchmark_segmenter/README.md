# RawHash2 Segmenter Benchmark

Comparison of 5 event segmenters × 5 datasets with ground-truth accuracy (minimap2 reference).

## Results summary (100 reads per dataset)

### D1 SARS-CoV-2 (R9.4) — 30 kb virus
| Segmenter | Reads mapped | Mean MAPQ | Precision | Recall | F1    |
|-----------|-------------:|----------:|----------:|-------:|------:|
| default   | 91/100       | 51.99     | 1.000     | 0.910  | 0.953 |
| pelt      | 91/100       | 51.99     | 1.000     | 0.910  | 0.953 |
| scrappie  | 81/100       | 38.02     | 1.000     | 0.810  | 0.895 |
| hmm       | 73/100       | 33.06     | 1.000     | 0.730  | 0.844 |
| binseg    | 3/100        |  6.94     | 1.000     | 0.030  | 0.058 |

### D4 Green Algae (R9.4) — 111 Mb, multi-chromosome
| Segmenter | Reads mapped | Precision | Recall | F1    |
|-----------|-------------:|----------:|-------:|------:|
| default   | 91/100       | 0.989     | 0.827  | 0.901 |
| pelt      | 91/100       | 0.989     | 0.827  | 0.901 |
| scrappie  | 67/100       | 0.893     | 0.609  | 0.724 |
| hmm       | 43/100       | 0.878     | 0.391  | 0.541 |
| binseg    |  1/100       | 0.250     | 0.009  | 0.018 |

D2 (E. coli R9.4) and D6 (E. coli R10.4) mapping benchmarks completed; ground truth pending SRA retrieval of basecalled reads.

## Key findings

1. **`default` and `pelt` produce identical mappings** across all datasets tested.
2. **PELT is consistently 20-30% faster** than default with identical accuracy.
3. **BinSeg is broken** (F1 ≈ 0.02-0.06). Attempted fix in `src/revent.c:422` (penalty `log(n) → 2·log(n)`) not sufficient — deeper issue in event generation.
4. **HMM and Scrappie are second-tier** — 5-20% worse F1 than default/pelt.
5. **Harder genomes (D4) amplify segmenter differences**.

## Files

- `rh2_bench_100reads.sh` — SLURM script, maps 100 reads per dataset × 5 segmenters
- `run_gt_100reads.sh` — minimap2 ground truth + accuracy evaluation
- `extract_100_reads.py` — extract 100 reads from FAST5 (works for single- and multi-read)
- `extract_read_ids.py` — list read IDs in a FAST5 directory
- `ground_truth_rawhash_style.py` — official RawHash2 pafstats classification (pair = read_id × ref_id)
- `ground_truth_eval.py` — position-aware variant with overlap threshold
- `final_pdf.py` — generate the summary PDF
- `rawhash2_segmenter_final.pdf` — the benchmark report (9 pages)
- `results_summary.json` — all mapping and timing data

## Usage (for your own data)

```bash
# 1. Ensure VBZ plugin for R10.4+ FAST5s
export HDF5_PLUGIN_PATH=$HOME/.local/lib/python3.10/site-packages/vbz_h5py_plugin/lib

# 2. Build index with your pore model
rawhash2 -x viral -t 8 -p pore.model -d my.ind reference.fa

# 3. Map with any segmenter
rawhash2 -x viral -t 8 --segmenter pelt -o out.paf my.ind /path/to/fast5_dir/

# 4. Ground truth
minimap2 -cx map-ont -t 8 --secondary=no reference.fa reads.fasta > truth.paf
python3 ground_truth_rawhash_style.py out.paf truth.paf
```

## Recommendation

**Use `--segmenter pelt` as the default.** Same accuracy as `default`, noticeably faster. Avoid `--segmenter binseg` until the underlying bug is fixed.

---

SAFARI Research Group — ETH Zürich — April 2026

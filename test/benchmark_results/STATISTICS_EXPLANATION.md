# Statistics Explanation — RawHash2 Segmenter Benchmark

> How every number in `rawhash2_segmenter_presentation.pdf` is computed.

---

## 1. Inputs used by the benchmark

| Input | Source | Purpose |
|---|---|---|
| `ref.fa` | Reference genome (FASTA) | Index target |
| `fast5_files/*.fast5` | Nanopore raw current recordings | Query reads |
| Pore model `*.model` | Oxford Nanopore k-mer expected currents | Signal simulation for indexing |
| `truth.paf` | minimap2 against basecalled reads | Ground-truth "correct" mapping positions |

All 5 segmenters receive **exactly the same** inputs, presets, thread count, and RawHash2 binary — the only knob changed is `--segmenter {default,hmm,pelt,binseg,scrappie}`.

---

## 2. Mapping Output Statistics (per segmenter, per dataset)

Extracted by parsing each `mappings.paf` file line-by-line. PAF columns used:

| Col | Name | Used for |
|---|---|---|
| 1 | Query name (read ID) | `unique_reads` |
| 10 | Residue matches | `mean_residue_matches` |
| 11 | Alignment block length | `mean_alignment_block_len` |
| 12 | MAPQ | `mapq_ge_X`, `mean_mapq`, `median_mapq` |

### `total_mappings`
> Total number of non-header lines in `mappings.paf`.

```python
total = sum(1 for line in open(paf))
```

### `unique_reads`
> Number of distinct read-IDs that produced at least one mapping.
> A read may appear multiple times (different chunks / secondary chains).

```python
unique = len({line.split('\t')[0] for line in open(paf)})
```

### `mapq_ge_X` (X in {1, 10, 30, 60})
> Count of mappings whose MAPQ is at least X.
> MAPQ ranges 0..60; higher = RawHash2 is more confident about the mapping location.

```python
mapq_ge_X = sum(1 for row in rows if int(row[11]) >= X)
```

### `mean_mapq`, `median_mapq`
> Arithmetic mean / median of MAPQ values across all mappings.

### `mean_alignment_block_len`
> Mean of column 11 across all mappings.
> The alignment block length is the number of events covered by the chain
> (including gaps) — essentially "how long was the match".

### `mean_residue_matches`
> Mean of column 10 — number of anchors in the chain that "matched"
> between read and reference.

---

## 3. Ground-Truth Accuracy (Precision / Recall / F1)

This is the most important block because it tells us if a mapping is **actually correct**, not just "confident".

### Setup

1. Basecall all FAST5 reads with Dorado → FASTA.
2. Map basecalled FASTA to `ref.fa` with `minimap2 -ax map-ont` → `truth.paf`.
3. For each read ID, `truth.paf` gives the canonical "ground-truth" region `(ref_start, ref_end, strand)`.
4. Compare each RawHash2 mapping against the ground truth by **overlap**.

### Overlap criterion

A RawHash2 mapping `(r_start_rh, r_end_rh)` is counted as **correct** (TP) if it
overlaps the ground-truth interval `(r_start_gt, r_end_gt)` by at least
`overlap_threshold = 10%` of the union length AND the strand matches.

```python
overlap = max(0, min(r_end_rh, r_end_gt) - max(r_start_rh, r_start_gt))
union  = max(r_end_rh, r_end_gt) - min(r_start_rh, r_start_gt)
correct = (overlap / union) >= 0.10 and strand_match
```

### Counts

| Metric | Definition |
|---|---|
| `n_truth_reads` | Number of reads that appear in `truth.paf` (minimap2 mapped them) |
| `n_query_reads` | Number of reads that RawHash2 mapped |
| `n_both_mapped` | Intersection of the two sets |
| **TP** (True Positive) | Reads where RawHash2 mapping matches truth |
| **FP** (False Positive) | Reads mapped by RawHash2 but to wrong location |
| **FN** (False Negative) | Reads mapped by truth but NOT by RawHash2 |
| `correct_position` | Same as TP |
| `wrong_position` | RawHash2 mapped, but location disagrees with truth (= FP) |

### Derived metrics

```
Precision = TP / (TP + FP)       "of the reads RawHash2 mapped, how many were correct?"
Recall    = TP / (TP + FN)       "of the mappable reads, how many did RawHash2 find?"
F1        = 2 * P * R / (P + R)  harmonic mean — one number that balances both
```

### MAPQ-split diagnostics

| Metric | Definition |
|---|---|
| `mean_mapq_correct` | Mean MAPQ of TP mappings. High value = RawHash2 is confident when it's right. |
| `mean_mapq_wrong` | Mean MAPQ of FP mappings. Low value = RawHash2 "knows" when it's unsure. |

A good segmenter has a wide gap between these two — the MAPQ signal is then a useful filter.

---

## 4. Performance Statistics

Collected by wrapping each RawHash2 call with `/usr/bin/time -v`. That produces a
standard stderr block like:

```
Elapsed (wall clock) time (h:mm:ss or m:ss): 13:33.41
User time (seconds): 1289.95
System time (seconds): 106.63
Maximum resident set size (kbytes): 4793246
```

### Parsed fields

| JSON key | Regex pattern | What it measures |
|---|---|---|
| `wall_clock_sec` | `Elapsed \(wall clock\) time.*: (.+)` | Real elapsed seconds (start → end) |
| `user_time_sec` | `User time \(seconds\): ([\d.]+)` | Summed CPU-user time over all threads |
| `sys_time_sec`  | `System time \(seconds\): ([\d.]+)` | Kernel-side CPU time |
| `max_rss_mb`    | `Maximum resident set size \(kbytes\): (\d+)` | Peak resident memory, kB → MB |

### Derived

```
throughput_reads_per_sec = n_reads / wall_clock_sec
cpu_efficiency           = user_time_sec / (wall_clock_sec * threads)  [0..1]
```

We run all segmenters with **8 threads** (`-t 8`) — same for every run — so
throughput comparisons are fair.

### Index-build vs mapping

We collect two separate `/usr/bin/time -v` runs per (dataset, segmenter):

1. **Indexing**: `rawhash2 -d index.ind ref.fa` — build the hash-table
2. **Mapping**: `rawhash2 -o mappings.paf index.ind subset/` — map all reads

In the PDF we report both as "Index Time (s)" and "Map Wall Time (s)".
Indexing typically takes < 1 s; mapping dominates.

### `index_size_bytes`
> `stat --format=%s index.ind` — size of the binary index on disk.
> Reported in the PDF as KB.

---

## 5. Cross-Dataset Aggregation

The final two tables average per-segmenter metrics across the 4 datasets
(D1, D2, D3, D6) with equal weight:

```
Mean_F1       = (F1_D1 + F1_D2 + F1_D3 + F1_D6) / 4
Mean_MAPQ     = ...
Mean_Tput     = ...
```

The "Final Ranking" table sorts segmenters by `Mean F1` — that is our primary
decision metric because it captures **both** whether reads are mapped (recall)
**and** whether they are mapped correctly (precision).

---

## 6. What each statistic tells you at a glance

| Metric | Direction | Purpose |
|---|---|---|
| Total mappings | informational | Volume of output |
| MAPQ >= 60 count | higher is better | How many "confident" mappings |
| Mean MAPQ | higher is better | Average confidence |
| **Precision** | higher is better | Are mappings correct? |
| **Recall** | higher is better | Are we finding the mappings? |
| **F1** | higher is better | Balanced overall accuracy |
| Wall clock time | lower is better | Real-world latency |
| Peak RSS | lower is better | Memory footprint |
| Throughput | higher is better | Scalability proxy |
| Index size | lower is better | Storage cost |

**Rule of thumb for presentation:** lead with F1, then throughput, then memory.

---

## 7. Why ground-truth accuracy differs from MAPQ

- MAPQ is RawHash2's **self-assessment**: "how confident am I?"
- F1 is the **external truth check**: "was I actually right?"

A segmenter can have high MAPQ but low F1 (overconfident) or low MAPQ but high
F1 (cautious but correct). Looking at both tells you whether the segmenter's
confidence signal is calibrated.

On D1 the Default segmenter has mean_mapq_correct = 52.5 vs mean_mapq_wrong = 9.9.
That gap of ~43 MAPQ units means a simple `MAPQ >= 10` filter removes most FPs
while keeping most TPs — the segmenter is producing a *useful* confidence signal.

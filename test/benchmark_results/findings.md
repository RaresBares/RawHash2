# Final Findings — Nanopore Segmenter Benchmark

## Setup
- 10 segmenters across 5 algorithmic families
- 7 datasets (D7 excluded — FAIL reads invalidate truth):
  - R9.4: D1 SARS-CoV-2, D2 E. coli, D3 yeast, D4 algae
  - R10.4: D6 E. coli
  - R10.4.1: RB_ecoli, RB_dmel, RB_hsap (full GRCh38)
- 3-fold cross-validation for HMM (12 configs) + BOCD (6 configs) per dataset
- Train-test gap = 0.000 across all CV folds (no overfitting)

## Per-dataset best F1 (segmenter / config):

| Dataset | Best F1 | Best segmenter | Default F1 | Δ |
|---|---|---|---|---|
| D1 | 0.968 | default | 0.968 | +0.000 |
| D2 | 0.530 | hmm | 0.410 | +0.120 |
| D3 | 0.545 | hmm | 0.295 | +0.250 |
| D4 | 0.373 | mad | 0.237 | +0.137 |
| D6 | 0.496 | binseg | 0.484 | +0.012 |
| RB_ecoli | 0.456 | hmm | 0.419 | +0.037 |
| RB_dmel | 0.482 | window | 0.480 | +0.002 |

## Per-segmenter sweet spots (where each is best):

- **hmm**: D2 (F1=0.53), D3 (F1=0.55), RB_ecoli (F1=0.46)
- **default**: D1 (F1=0.97)
- **binseg**: D6 (F1=0.50)
- **window**: RB_dmel (F1=0.48)
- **mad**: D4 (F1=0.37)

## Speed champions (≥95% F1 of best, fastest):

| Dataset | Default wall | Fast segmenter (≥95% F1) | wall | Speedup |
|---|---|---|---|---|
| D1 | 1359s | window (F1=0.953) | 651s | 2.09× |
| D2 | 549s | hmm (F1=0.530) | 439s | 1.25× |
| D3 | 6s | hmm (F1=0.545) | 4s | 1.52× |
| D4 | 13s | mad (F1=0.373) | 16s | 0.82× |
| D6 | 1332s | binseg (F1=0.496) | 1203s | 1.11× |
| RB_ecoli | 55s | hmm (F1=0.456) | 76s | 0.71× |
| RB_dmel | 366s | window (F1=0.482) | 339s | 1.08× |

## CV-tuned best config per (dataset, segmenter)

- **D1/hmm**: RH2_HMM_STATES=6 RH2_HMM_STAY=0.96 → F1_test = 0.8675 ± 0.0004
- **D1/bocd**: RH2_BOCD_EXP_LEN=600 → F1_test = 0.5068 ± 0.0010
- **D2/bocd**: RH2_BOCD_EXP_LEN=600 → F1_test = 0.3877 ± 0.0004
- **D2/hmm**: RH2_HMM_STATES=4 RH2_HMM_STAY=0.93 → F1_test = 0.5381 ± 0.0015
- **D4/bocd**: RH2_BOCD_EXP_LEN=600 → F1_test = 0.0034 ± 0.0002
- **D4/hmm**: RH2_HMM_STATES=5 RH2_HMM_STAY=0.96 → F1_test = 0.0058 ± 0.0008
- **D6/hmm**: RH2_HMM_STATES=3 RH2_HMM_STAY=0.93 → F1_test = 0.0162 ± 0.0001
- **D6/bocd**: RH2_BOCD_EXP_LEN=600 → F1_test = 0.3187 ± 0.0009
- **RB_dmel/hmm**: RH2_HMM_STATES=6 RH2_HMM_STAY=0.96 → F1_test = 0.4286 ± 0.0064
- **RB_dmel/bocd**: RH2_BOCD_EXP_LEN=600 → F1_test = 0.3560 ± 0.0105
- **RB_ecoli/bocd**: RH2_BOCD_EXP_LEN=600 → F1_test = 0.3566 ± 0.0061
- **RB_ecoli/hmm**: RH2_HMM_STATES=6 RH2_HMM_STAY=0.96 → F1_test = 0.5040 ± 0.0036

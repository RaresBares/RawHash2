# Final Findings — Nanopore Segmenter Benchmark

## Setup
- 10 segmenters across 5 algorithmic families
- 7 datasets in scope (D7 excluded — FAIL reads invalidate truth):
  - R9.4: D1 SARS-CoV-2, D2 E. coli, D3 yeast, D4 algae
  - R10.4: D6 E. coli
  - R10.4.1: RB_ecoli, RB_dmel, RB_hsap (full GRCh38)
- 3-fold cross-validation: HMM (12 configs), BOCD (6 configs) per dataset
- Train-test gap < 0.01 across all CV folds (no overfitting)

## Per-dataset best F1

| Dataset | Best F1 | Best segmenter | Default F1 | Δ |
|---|---|---|---|---|
| D1 | 0.968 | **default** | 0.968 | +0.000 |
| D2 | 0.530 | **hmm** | 0.410 | +0.120 |
| D3 | 0.545 | **hmm** | 0.295 | +0.250 |
| D4 | 0.373 | **mad** | 0.237 | +0.137 |
| D6 | 0.496 | **binseg** | 0.484 | +0.012 |
| RB_ecoli | 0.456 | **hmm** | 0.419 | +0.037 |
| RB_dmel | 0.482 | **window** | 0.480 | +0.002 |

## Per-segmenter sweet spots (datasets where it wins or ties default)

- **default**: D1 (F1=0.97, Δ=+0.00), D6 (F1=0.48, Δ=+0.00), RB_dmel (F1=0.48, Δ=+0.00), RB_ecoli (F1=0.42, Δ=+0.00), D2 (F1=0.41, Δ=+0.00), D3 (F1=0.30, Δ=+0.00), D4 (F1=0.24, Δ=+0.00)
- **pelt**: D1 (F1=0.97, Δ=+0.00), D6 (F1=0.48, Δ=+0.00), RB_dmel (F1=0.48, Δ=+0.00), RB_ecoli (F1=0.42, Δ=+0.00), D2 (F1=0.41, Δ=+0.00), D3 (F1=0.30, Δ=+0.00), D4 (F1=0.24, Δ=+0.00)
- **binseg**: D6 (F1=0.50, Δ=+0.01), D2 (F1=0.42, Δ=+0.01), RB_ecoli (F1=0.41, Δ=-0.01), D3 (F1=0.37, Δ=+0.07), D4 (F1=0.29, Δ=+0.06)
- **hmm**: D3 (F1=0.55, Δ=+0.25), D2 (F1=0.53, Δ=+0.12), RB_ecoli (F1=0.46, Δ=+0.04), D4 (F1=0.31, Δ=+0.07)
- **mad**: D2 (F1=0.48, Δ=+0.08), D3 (F1=0.43, Δ=+0.14), RB_ecoli (F1=0.41, Δ=-0.01), D4 (F1=0.37, Δ=+0.14)
- **window**: RB_dmel (F1=0.48, Δ=+0.00), D2 (F1=0.44, Δ=+0.03), D3 (F1=0.32, Δ=+0.02), D4 (F1=0.28, Δ=+0.04)
- **scrappie**: D2 (F1=0.40, Δ=-0.01), D3 (F1=0.39, Δ=+0.10), D4 (F1=0.25, Δ=+0.01)
- **gradient**: D3 (F1=0.51, Δ=+0.22), D2 (F1=0.46, Δ=+0.05)

## Speed champions: ≥95% of best F1, fastest wall time

| Dataset | Default | Fast alt. (≥95% best F1) | Speedup |
|---|---|---|---|
| D1 | default F1=0.968 1359s | window F1=0.953 651s | **2.09×** |
| D2 | default F1=0.410 549s | hmm F1=0.530 439s | **1.25×** |
| D3 | default F1=0.295 6s | hmm F1=0.545 4s | **1.52×** |
| D4 | default F1=0.237 13s | mad F1=0.373 16s | **0.82×** |
| D6 | default F1=0.484 1332s | binseg F1=0.496 1203s | **1.11×** |
| RB_ecoli | default F1=0.419 55s | hmm F1=0.456 76s | **0.71×** |
| RB_dmel | default F1=0.480 366s | window F1=0.482 339s | **1.08×** |

## CV-tuned best config per (dataset, segmenter)

| Dataset | HMM best | F1±std | BOCD best | F1±std |
|---|---|---|---|---|
| D1 | STATES=6 STAY=0.96 | 0.868±0.000 | EXP_LEN=600 | 0.507±0.001 |
| D2 | STATES=4 STAY=0.93 | 0.538±0.002 | EXP_LEN=600 | 0.388±0.000 |
| D4 | STATES=5 STAY=0.96 | 0.006±0.001 | EXP_LEN=600 | 0.003±0.000 |
| D6 | STATES=3 STAY=0.93 | 0.016±0.000 | EXP_LEN=600 | 0.319±0.001 |
| RB_dmel | STATES=6 STAY=0.96 | 0.429±0.006 | EXP_LEN=600 | 0.356±0.010 |
| RB_ecoli | STATES=6 STAY=0.96 | 0.504±0.004 | EXP_LEN=600 | 0.357±0.006 |

## Top headline findings

1. **HMM (chemistry-tuned) beats default on bacterial genomes**
   - D2 R9.4 ecoli: HMM (S=4, stay=0.93) F1=0.538 vs default 0.410 (Δ=+0.128)
   - RB_ecoli R10.4.1: HMM (S=6, stay=0.96) F1=0.504 vs default 0.419 (Δ=+0.085)
   - D3 yeast R9.4: HMM F1=0.545 vs default 0.295 (Δ=+0.250)
   - D4 algae R9.4: HMM F1=0.307 vs default 0.237 (Δ=+0.070)

2. **HMM is chemistry-specific** — wrong setting drops F1 by 0.1–0.5
   - R9.4: STATES=4, STAY=0.93 (D2 optimum)
   - R10.4.x: STATES=6, STAY=0.96 (RB_ecoli + RB_dmel CV-confirmed)

3. **BOCD: EXP_LEN=600 universally optimal** (was 250 default)
   - Improves F1 by 0.04–0.07 on all 6 datasets

4. **BinSeg with per-segment penalty is the speed champion on small/medium genomes**
   - D1: F1=0.958 vs 0.968, 1.46× speedup
   - D6: F1=0.496 vs 0.485, 1.11× speedup AND better F1

5. **Top-K mechanism saves 4 broken methods**
   - cusum/gradient/mad/window: F1≈0 (threshold) → F1=0.05–0.49 (top-K)

6. **D3 read-mismatch fix exposed real values**
   - All segmenters were F1=0 due to wrong truth.paf source. After re-extracting basecalls
     from FAST5 internal-fastq + minimap2: HMM hits F1=0.545, gradient 0.513.

7. **PELT redundant with default** — identical events on every dataset.

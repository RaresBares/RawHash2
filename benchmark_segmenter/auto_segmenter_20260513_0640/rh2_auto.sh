#!/bin/bash
# rh2_auto.sh - run rawhash2 with the rule-picked segmenter for one dataset.
#
# Mirrors the per-segmenter call in bench_v4_full.sh so wall times stay
# comparable to v4_aggregated.json. Sets the per-segmenter tuned env vars
# (HMM states/stay, BOCD exp_len, PELT pen/min_size) before launch, and
# enables RH2_TOPK_AUTO when the rule picks an "auto" mode.
#
# Required env: DS, REF, F5, PORE, OUTDIR, CHEM (R9.4|R10.4|R10.4.1), EXTRA
# Optional:     T (threads, default 8)

set -euo pipefail
DS=${DS:?}; REF=${REF:?}; F5=${F5:?}; PORE=${PORE:?}; OUTDIR=${OUTDIR:?}
CHEM=${CHEM:?}; EXTRA=${EXTRA:-}
T=${T:-8}
RH=${RH:-/home/rsahleanu/rawhash2/bin/rawhash2}
PICK=${PICK:-$(dirname "$0")/pick_segmenter.py}

mkdir -p "$OUTDIR"

# Ask the rule for (segmenter, mode); writes shell-eval'able exports.
eval "$(python3 "$PICK" --ref "$REF" --chemistry "$CHEM" --shell)"
echo "[$DS] rule -> segmenter=$AUTO_SEG mode=$AUTO_MODE (ref_mb=$AUTO_REF_MB gc=$AUTO_REF_GC)"

# Per-segmenter tuned configs (same as bench_v4_full.sh).
case "$AUTO_SEG" in
    hmm)
        case "$DS" in
            D1)         export RH2_HMM_STATES=6; export RH2_HMM_STAY=0.96 ;;
            D2)         export RH2_HMM_STATES=4; export RH2_HMM_STAY=0.93 ;;
            D3)         export RH2_HMM_STATES=6; export RH2_HMM_STAY=0.93 ;;
            D4)         export RH2_HMM_STATES=5; export RH2_HMM_STAY=0.96 ;;
            *)          export RH2_HMM_STATES=6; export RH2_HMM_STAY=0.96 ;;
        esac
        ;;
    bocd)
        export RH2_BOCD_EXP_LEN=600
        ;;
    pelt|pelt_cuda)
        export RH2_PELT_PEN_MULT=0.05
        export RH2_PELT_MIN_SIZE=3
        ;;
esac

IDX="$OUTDIR/idx_${AUTO_SEG}_${AUTO_MODE}.ind"
PAF="$OUTDIR/auto_${AUTO_SEG}_${AUTO_MODE}.paf"

if [ ! -s "$IDX" ]; then
    echo "[$DS] index..."
    $RH -x sensitive -t $T -p "$PORE" $EXTRA $AUTO_SEG_FLAG \
        -d "$IDX" "$REF" 2>>"$OUTDIR/idx.log"
fi

echo "[$DS] map..."
T0=$(date +%s.%N)
$RH -x sensitive -t $T $EXTRA $AUTO_SEG_FLAG -o "$PAF" "$IDX" "$F5" \
    2>>"$OUTDIR/map.log" || true
T1=$(date +%s.%N)
WALL=$(echo "$T1-$T0" | bc)

echo "[$DS] done seg=$AUTO_SEG mode=$AUTO_MODE wall=${WALL}s -> $PAF"
echo "$DS,$AUTO_SEG,$AUTO_MODE,$WALL" >> "$OUTDIR/auto_results.csv"

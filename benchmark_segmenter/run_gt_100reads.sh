#!/bin/bash
#SBATCH --job-name=rh2_gt100
#SBATCH --output=/mnt/galactica/rsahleanu/seq_benchmark/results/rh2_gt100.out
#SBATCH --error=/mnt/galactica/rsahleanu/seq_benchmark/results/rh2_gt100.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu_part
set -euo pipefail

source /home/rsahleanu/miniconda3/etc/profile.d/conda.sh
conda activate nanopore
export HDF5_PLUGIN_PATH=/home/rsahleanu/.local/lib/python3.10/site-packages/vbz_h5py_plugin/lib

W=/mnt/galactica/rsahleanu/seq_benchmark
R=$W/results/rh2_100reads
GT_SCRIPT=$W/scripts/ground_truth_eval.py
ID_SCRIPT=$W/scripts/extract_read_ids.py
SEGS=(default hmm pelt binseg scrappie)

echo "=== GT Eval on 100-read subsets - $(date) ==="

for DS_DIR in d1_sars_r94 d2_ecoli_r94 d3_yeast_r94 d4_green_algae_r94 d6_ecoli_r104; do
    case $DS_DIR in
        d1_sars_r94)       DS=D1;;
        d2_ecoli_r94)      DS=D2;;
        d3_yeast_r94)      DS=D3;;
        d4_green_algae_r94) DS=D4;;
        d6_ecoli_r104)     DS=D6;;
    esac

    READS=$W/datasets/$DS_DIR/reads.fasta
    REF=$W/datasets/$DS_DIR/ref.fa
    SUBSET_F5=$W/datasets/$DS_DIR/fast5_100reads

    if [ ! -s "$READS" ]; then echo "  [SKIP] $DS: no reads.fasta"; continue; fi
    if [ ! -s "$REF" ]; then echo "  [SKIP] $DS: no ref.fa"; continue; fi
    if [ ! -d "$R/$DS" ]; then echo "  [SKIP] $DS: no benchmark results"; continue; fi
    if [ ! -d "$SUBSET_F5" ]; then echo "  [SKIP] $DS: no fast5_100reads"; continue; fi

    echo ""
    echo "=== $DS ==="

    # Extract read IDs of the 100-read subset
    ID_LIST=$R/$DS/subset_read_ids.txt
    python3 "$ID_SCRIPT" "$SUBSET_F5" "$ID_LIST" 2>&1
    NIDS=$(wc -l < "$ID_LIST")
    echo "  Subset read IDs: $NIDS"

    # Generate truth.paf if needed (full reads.fasta vs ref)
    GT=$R/$DS/truth.paf
    if [ ! -s "$GT" ]; then
        echo "  [$DS] minimap2 ..."
        minimap2 -cx map-ont -t 8 --secondary=no "$REF" "$READS" 2>/dev/null > "$GT" || true
    fi
    NGT=$(wc -l < "$GT")
    echo "  truth.paf: $NGT mappings"

    for S in "${SEGS[@]}"; do
        QUERY=$R/$DS/$S/mappings.paf
        OUT=$R/$DS/$S/accuracy.json
        if [ -s "$QUERY" ]; then
            python3 "$GT_SCRIPT" "$GT" "$QUERY" 0.1 "$ID_LIST" > "$OUT" 2>/dev/null || echo '{}' > "$OUT"
            P=$(python3 -c "import json;print(json.load(open('$OUT')).get('precision',0))")
            REC=$(python3 -c "import json;print(json.load(open('$OUT')).get('recall',0))")
            F1=$(python3 -c "import json;print(json.load(open('$OUT')).get('f1',0))")
            echo "  [$DS/$S] P=$P R=$REC F1=$F1"
        fi
    done
done

echo ""
echo "=== Regenerate final PDF with accuracy ==="
python3 $W/scripts/rh2_collect_and_plot_v3.py $R $R/results_with_gt.json

echo "=== DONE $(date) ==="
echo "PDF: $R/rawhash2_segmenter_comparison_v3.pdf"

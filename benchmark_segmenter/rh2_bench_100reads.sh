#!/bin/bash
#SBATCH --job-name=rh2_100
#SBATCH --output=/mnt/galactica/rsahleanu/seq_benchmark/results/rh2_100.out
#SBATCH --error=/mnt/galactica/rsahleanu/seq_benchmark/results/rh2_100.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu_part
set -euo pipefail

export HDF5_PLUGIN_PATH=/home/rsahleanu/.local/lib/python3.10/site-packages/vbz_h5py_plugin/lib

W=/mnt/galactica/rsahleanu/seq_benchmark
RH=$HOME/rawhash2/bin/rawhash2
R=$W/results/rh2_100reads
P94=$HOME/rawhash2/extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model
P104=$HOME/rawhash2/extern/kmer_models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt
T=8
SEGS=(default hmm pelt binseg scrappie)

mkdir -p $R
echo "=== RawHash2 100-Read Benchmark - $(date) ==="

# Format: DS:dirname:preset:pore_model:extra_flags
CONFIGS=(
    "D1:d1_sars_r94:viral:$P94:"
    "D2:d2_ecoli_r94:sensitive:$P94:"
    "D3:d3_yeast_r94:sensitive:$P94:"
    "D4:d4_green_algae_r94:sensitive:$P94:"
    "D6:d6_ecoli_r104:sensitive:$P104:--r10"
)

for entry in "${CONFIGS[@]}"; do
    IFS=':' read -r DS DIRNAME PX PM EX <<< "$entry"
    DD=$W/datasets/$DIRNAME
    SUB=$DD/fast5_100reads

    # Skip if extraction not done
    if [ ! -d "$SUB" ] || [ -z "$(ls $SUB/*.fast5 2>/dev/null)" ]; then
        echo "=== [SKIP] $DS: no fast5_100reads at $SUB ==="
        continue
    fi

    NA=$(ls $SUB/*.fast5 | wc -l)
    echo ""
    echo "======================================================================"
    echo "  $DS (preset=$PX $EX) — $NA FAST5 file(s) with ~100 reads"
    echo "======================================================================"

    rm -rf $R/$DS
    mkdir -p $R/$DS

    for S in "${SEGS[@]}"; do
        SD=$R/$DS/$S; mkdir -p $SD
        echo "  [$DS/$S] index..."
        /usr/bin/time -v $RH -x $PX $EX -t $T --segmenter $S -p "$PM" -d "$SD/index.ind" $DD/ref.fa 2>"$SD/index_time.txt" || true
        echo "  [$DS/$S] map..."
        /usr/bin/time -v $RH -x $PX $EX -t $T --segmenter $S -o "$SD/mappings.paf" "$SD/index.ind" "$SUB" 2>"$SD/map_time.txt" || true
        [ -f "$SD/mappings.paf" ] && echo "  [$DS/$S] $(wc -l<$SD/mappings.paf) maps" || echo "  [$DS/$S] FAIL"
    done

    # Link n_reads info (the extractor directory name says 100)
    mkdir -p $R/$DS/subset
    for f in $SUB/*.fast5; do ln -sf "$f" $R/$DS/subset/; done
done

echo ""
echo "=== Regenerate PDF ==="
python3 $W/scripts/rh2_collect_and_plot_v2.py $R $R/results.json

echo "=== DONE $(date) ==="

#!/bin/bash
#SBATCH --job-name=rh2_download
#SBATCH --partition=cpu_part
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/home/rsahleanu/rawhash2/test/benchmark_results/download_%j.out
#SBATCH --error=/home/rsahleanu/rawhash2/test/benchmark_results/download_%j.err

set -e
DATA_DIR="/home/rsahleanu/rawhash2/test/data"
echo "Download started: $(date)"

# ---- D2: E. coli R9.4 ----
echo ""
echo ">>> D2: E. coli R9.4"
if [ "$(ls "$DATA_DIR/d2_ecoli_r94/fast5_files/"*.fast5 2>/dev/null | wc -l)" -lt 50 ]; then
    mkdir -p "$DATA_DIR/d2_ecoli_r94/fast5_files"
    cd "$DATA_DIR/d2_ecoli_r94"
    echo "  Downloading and extracting FAST5 files..."
    wget -q -O- https://sra-pub-src-2.s3.amazonaws.com/ERR9127551/ecoli_r9.tar.gz.1 | tar xz || true
    find . -type f -name '*.fast5' -not -path './fast5_files/*' | head -200 | xargs -I{} mv {} fast5_files/ 2>/dev/null
    # Clean up extracted dirs
    find . -mindepth 1 -maxdepth 1 -type d -not -name fast5_files -exec rm -rf {} + 2>/dev/null
    echo "  D2 fast5 count: $(ls fast5_files/*.fast5 2>/dev/null | wc -l)"
else
    echo "  D2 already has enough files: $(ls "$DATA_DIR/d2_ecoli_r94/fast5_files/"*.fast5 2>/dev/null | wc -l)"
fi

# ---- D6: E. coli R10.4 ----
echo ""
echo ">>> D6: E. coli R10.4"
if [ "$(ls "$DATA_DIR/d6_ecoli_r104/fast5_files/"*.fast5 2>/dev/null | wc -l)" -lt 50 ]; then
    mkdir -p "$DATA_DIR/d6_ecoli_r104/fast5_files"
    cd "$DATA_DIR/d6_ecoli_r104"
    echo "  Downloading and extracting FAST5 files..."
    wget -q -O- https://sra-pub-src-2.s3.amazonaws.com/ERR9127552/Ecoli_r10.4.tar.gz.1 | tar xz || true
    find . -type f -name '*.fast5' -not -path './fast5_files/*' | head -200 | xargs -I{} mv {} fast5_files/ 2>/dev/null
    find . -mindepth 1 -maxdepth 1 -type d -not -name fast5_files -exec rm -rf {} + 2>/dev/null
    echo "  D6 fast5 count: $(ls fast5_files/*.fast5 2>/dev/null | wc -l)"
else
    echo "  D6 already has enough files: $(ls "$DATA_DIR/d6_ecoli_r104/fast5_files/"*.fast5 2>/dev/null | wc -l)"
fi

echo ""
echo "=== Download Summary ==="
for d in d2_ecoli_r94 d6_ecoli_r104; do
    n=$(ls "$DATA_DIR/$d/fast5_files/"*.fast5 2>/dev/null | wc -l)
    echo "  $d: $n fast5 files"
done
echo "Download finished: $(date)"

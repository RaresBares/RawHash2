#!/bin/bash
set -e

# Download 6 selected datasets for segmenter benchmarking
# Only download FAST5 files + reference genomes (skip reads.fasta that needs external tools)
# We only need 50 reads per dataset

DATA_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DATA_DIR"

echo "=== Downloading 6 datasets for segmenter benchmarking ==="
echo "Start time: $(date)"

# ---- D1: SARS-CoV-2 R9.4 (small) ----
echo ""
echo ">>> D1: SARS-CoV-2 R9.4 (small genome ~30kb)"
if [ ! -d d1_sars-cov-2_r94/fast5_files ] || [ "$(ls d1_sars-cov-2_r94/fast5_files/*.fast5 2>/dev/null | wc -l)" -eq 0 ]; then
    mkdir -p d1_sars-cov-2_r94/fast5_files/
    cd d1_sars-cov-2_r94
    echo "  Downloading FAST5 files..."
    wget -q --show-progress -O- https://cadde.s3.climb.ac.uk/SP1-raw.tgz | tar -xz
    find ./SP1-fast5-mapped -type f -name '*.fast5' | xargs -I{} mv {} ./fast5_files/
    rm -rf SP1-fast5-mapped SP1-mapped.fastq README 2>/dev/null
    echo "  Downloading reference genome..."
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
    gunzip -f GCF_009858895.2_ASM985889v3_genomic.fna.gz
    mv GCF_009858895.2_ASM985889v3_genomic.fna ref.fa
    cd "$DATA_DIR"
    echo "  D1 done: $(ls d1_sars-cov-2_r94/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
else
    echo "  D1 already downloaded: $(ls d1_sars-cov-2_r94/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
fi

# ---- D7: S. aureus R10.4 (small) ----
echo ""
echo ">>> D7: S. aureus R10.4 (small genome ~2.9Mb)"
if [ ! -d d7_saureus_r104/fast5_files ] || [ "$(ls d7_saureus_r104/fast5_files/*.fast5 2>/dev/null | wc -l)" -eq 0 ]; then
    mkdir -p d7_saureus_r104/fast5_files/
    cd d7_saureus_r104/fast5_files
    echo "  Downloading FAST5 files..."
    wget -q --show-progress -O- https://sra-pub-src-1.s3.amazonaws.com/SRR21386013/S_aureus_JKD6159_ONT_R10.4_fast5.tar.gz.1 | tar xz
    # Move fast5 files from nested directories
    find . -type f -name '*.fast5' -not -path './*.fast5' -exec mv {} . \; 2>/dev/null
    find . -mindepth 1 -type d -exec rm -rf {} + 2>/dev/null
    cd "$DATA_DIR/d7_saureus_r104"
    echo "  Downloading reference genome..."
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/144/955/GCF_000144955.2_ASM14495v2/GCF_000144955.2_ASM14495v2_genomic.fna.gz
    gunzip -f GCF_000144955.2_ASM14495v2_genomic.fna.gz
    mv GCF_000144955.2_ASM14495v2_genomic.fna ref.fa
    cd "$DATA_DIR"
    echo "  D7 done: $(ls d7_saureus_r104/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
else
    echo "  D7 already downloaded: $(ls d7_saureus_r104/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
fi

# ---- D2: E. coli R9.4 (medium) ----
echo ""
echo ">>> D2: E. coli R9.4 (medium genome ~5Mb)"
if [ ! -d d2_ecoli_r94/fast5_files ] || [ "$(ls d2_ecoli_r94/fast5_files/*.fast5 2>/dev/null | wc -l)" -eq 0 ]; then
    mkdir -p d2_ecoli_r94/fast5_files/
    cd d2_ecoli_r94
    echo "  Downloading FAST5 files..."
    wget -q --show-progress -O- https://sra-pub-src-2.s3.amazonaws.com/ERR9127551/ecoli_r9.tar.gz.1 | tar xz
    mv r9/f5s/RefStrains210914_NK/f5s/barcode02/*.fast5 fast5_files/ 2>/dev/null
    rm -rf r9
    echo "  Downloading reference genome..."
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/445/GCA_000007445.1_ASM744v1/GCA_000007445.1_ASM744v1_genomic.fna.gz
    gunzip -f GCA_000007445.1_ASM744v1_genomic.fna.gz
    mv GCA_000007445.1_ASM744v1_genomic.fna ref.fa
    cd "$DATA_DIR"
    echo "  D2 done: $(ls d2_ecoli_r94/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
else
    echo "  D2 already downloaded: $(ls d2_ecoli_r94/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
fi

# ---- D6: E. coli R10.4 (medium) ----
echo ""
echo ">>> D6: E. coli R10.4 (medium genome ~5Mb)"
if [ ! -d d6_ecoli_r104/fast5_files ] || [ "$(ls d6_ecoli_r104/fast5_files/*.fast5 2>/dev/null | wc -l)" -eq 0 ]; then
    mkdir -p d6_ecoli_r104/fast5_files/
    cd d6_ecoli_r104
    echo "  Downloading FAST5 files..."
    wget -q --show-progress -O- https://sra-pub-src-2.s3.amazonaws.com/ERR9127552/Ecoli_r10.4.tar.gz.1 | tar xz
    find ./mnt -type f -name '*.fast5' -exec mv {} fast5_files/ \; 2>/dev/null
    rm -rf ./mnt
    echo "  Downloading reference genome..."
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/445/GCA_000007445.1_ASM744v1/GCA_000007445.1_ASM744v1_genomic.fna.gz
    gunzip -f GCA_000007445.1_ASM744v1_genomic.fna.gz
    mv GCA_000007445.1_ASM744v1_genomic.fna ref.fa
    cd "$DATA_DIR"
    echo "  D6 done: $(ls d6_ecoli_r104/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
else
    echo "  D6 already downloaded: $(ls d6_ecoli_r104/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
fi

# ---- D3: Yeast R9.4 (large) ----
echo ""
echo ">>> D3: Yeast R9.4 (large genome ~12Mb)"
if [ ! -d d3_yeast_r94/fast5_files ] || [ "$(ls d3_yeast_r94/fast5_files/*.fast5 2>/dev/null | wc -l)" -eq 0 ]; then
    mkdir -p d3_yeast_r94/fast5_files/
    cd d3_yeast_r94
    echo "  Downloading FAST5 files (this may take a while)..."
    wget -q --show-progress -O- https://sra-pub-src-1.s3.amazonaws.com/SRR8648503/GLU1II_basecalled_fast5_1.tar.gz.1 | tar -xz
    # Only take first 200 files (we only need 50 reads, but take extras for safety)
    find ./GLU1II_basecalled_fast5_1 -type f -name '*.fast5' | head -200 | xargs -I{} mv {} ./fast5_files/
    rm -rf GLU1II_basecalled_fast5_1
    echo "  Downloading reference genome..."
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
    gunzip -f sacCer3.fa.gz
    mv sacCer3.fa ref.fa
    cd "$DATA_DIR"
    echo "  D3 done: $(ls d3_yeast_r94/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
else
    echo "  D3 already downloaded: $(ls d3_yeast_r94/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
fi

# ---- D4: Green Algae R9.4 (large) ----
echo ""
echo ">>> D4: Green Algae R9.4 (large genome ~111Mb)"
if [ ! -d d4_green_algae_r94/fast5_files ] || [ "$(ls d4_green_algae_r94/fast5_files/*.fast5 2>/dev/null | wc -l)" -eq 0 ]; then
    mkdir -p d4_green_algae_r94/fast5_files/
    cd d4_green_algae_r94
    echo "  Downloading FAST5 files (this may take a while)..."
    wget -q --show-progress -O- https://sra-pub-src-2.s3.amazonaws.com/ERR3237140/Chlamydomonas_0.tar.gz.1 | tar xz
    find ./Chlamydomonas_0/reads/downloads/pass/ -type f -name '*.fast5' | head -200 | xargs -I{} mv {} fast5_files/
    rm -rf Chlamydomonas_0
    echo "  Downloading reference genome..."
    wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/595/GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5/GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.fna.gz
    gunzip -f GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.fna.gz
    mv GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.fna ref.fa
    cd "$DATA_DIR"
    echo "  D4 done: $(ls d4_green_algae_r94/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
else
    echo "  D4 already downloaded: $(ls d4_green_algae_r94/fast5_files/*.fast5 2>/dev/null | wc -l) fast5 files"
fi

echo ""
echo "=== Download Summary ==="
for d in d1_sars-cov-2_r94 d7_saureus_r104 d2_ecoli_r94 d6_ecoli_r104 d3_yeast_r94 d4_green_algae_r94; do
    n=$(ls "$DATA_DIR/$d/fast5_files/"*.fast5 2>/dev/null | wc -l)
    ref_size=$(stat --format=%s "$DATA_DIR/$d/ref.fa" 2>/dev/null || echo "0")
    echo "  $d: $n fast5 files, ref.fa = $(echo "scale=1; $ref_size / 1048576" | bc)MB"
done
echo ""
echo "End time: $(date)"
echo "=== All downloads complete ==="

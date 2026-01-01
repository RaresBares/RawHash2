#!/bin/bash

DATASET="d5_human_na12878_r94"
SIGNALS="../../../data/${DATASET}/pod5_files/"
READS="../../../data/${DATASET}/reads.fasta"
REF="../../../data/${DATASET}/ref.fa"
PRESET="ava"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"

OUTDIR="rawsamble"
THREAD=$1 #64

mkdir -p ${OUTDIR}

INDEX="${OUTDIR}/${DATASET}_rawsamble_${PRESET}.ind"

/usr/bin/time -vpo "${OUTDIR}/${DATASET}_rawsamble_index_${PRESET}.time" rawhash2 -x ${PRESET} -t ${THREAD} -p "${PORE}" -d ${INDEX} ${SIGNALS} > "${OUTDIR}/${DATASET}_rawsamble_index_${PRESET}.out" 2> "${OUTDIR}/${DATASET}_rawsamble_index_${PRESET}.err"

/usr/bin/time -vpo "${OUTDIR}/${DATASET}_rawsamble_map_${PRESET}.time" rawhash2 -x ${PRESET} -t ${THREAD} ${INDEX} ${SIGNALS} > "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.paf" 2> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.err"

miniasm "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.paf" > "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.gfa"

echo "GFA Stats from Rawsamble mappings:" >> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results"
bash ../../../scripts/analyze_gfa.sh ${OUTDIR}/${DATASET}_rawsamble_${PRESET}.gfa > "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results"
python ../../../scripts/compute_aun.py "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.gfa"

/usr/bin/time -vpo "${OUTDIR}/mm2_overlaps.time" minimap2 -t ${THREAD} -x ava-ont --for-only ${READS} ${READS} > "${OUTDIR}/mm2_overlaps.paf"
miniasm "${OUTDIR}/mm2_overlaps.paf" > "${OUTDIR}/mm2_overlaps.gfa"
echo "GFA Stats from minimap2 mappings:" >> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results"
bash ../../../scripts/analyze_gfa.sh  "${OUTDIR}/mm2_overlaps.gfa" >> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results"
python ../../../scripts/compute_aun.py "${OUTDIR}/mm2_overlaps.gfa"

python3 ../../../scripts/pafstats.py "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.paf" "${OUTDIR}/mm2_overlaps.paf" 2>> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results"

bash ../../../scripts/run_minimap2_multimap.sh "${OUTDIR}" ${READS} ${REF} ${THREAD}
echo "Rawsamble chained read percentage:" >> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results"
python ../../../scripts/evaluate_gfa.py "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.gfa" "${OUTDIR}/true_mappings.paf" >> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results"
echo "Minimap2 chained read percentage:" >> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results"
python ../../../scripts/evaluate_gfa.py "${OUTDIR}/mm2_overlaps.gfa" "${OUTDIR}/true_mappings.paf" >> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results" >> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.results"

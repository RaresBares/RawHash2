#!/bin/bash

DATASET="d9_ecoli_r1041"
SIGNALS="../../../data/${DATASET}/nonsplit_pod5_files/"
READS="../../../data/${DATASET}/reads_dorado.fasta"
REF="../../../data/${DATASET}/ref.fa"
PRESET="ava"
PORE="../../../../extern/local_kmer_models/uncalled_r1041_model_only_means.txt"

OUTDIR="rawsamble"
THREAD=$1 #64
PARAMS="--r10"

mkdir -p ${OUTDIR}

INDEX="${OUTDIR}/${DATASET}_rawsamble_${PRESET}.ind"

/usr/bin/time -vpo "${OUTDIR}/${DATASET}_rawsamble_index_${PRESET}.time" rawhash2_amd -x ${PRESET} -t ${THREAD} -p "${PORE}" -d ${INDEX} "${PARAMS}" ${SIGNALS} > "${OUTDIR}/${DATASET}_rawsamble_index_${PRESET}.out" 2> "${OUTDIR}/${DATASET}_rawsamble_index_${PRESET}.err"

/usr/bin/time -vpo "${OUTDIR}/${DATASET}_rawsamble_map_${PRESET}.time" rawhash2_amd -x ${PRESET} -t ${THREAD} "${PARAMS}" ${INDEX} ${SIGNALS} > "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.paf" 2> "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.err"

miniasm "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.paf" > "${OUTDIR}/${DATASET}_rawsamble_${PRESET}.gfa"

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

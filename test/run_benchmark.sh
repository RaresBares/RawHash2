#!/bin/bash
#SBATCH --job-name=rawhash2_segbench
#SBATCH --partition=cpu_part
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/rsahleanu/rawhash2/test/benchmark_results/slurm_%j.out
#SBATCH --error=/home/rsahleanu/rawhash2/test/benchmark_results/slurm_%j.err

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

mkdir -p /home/rsahleanu/rawhash2/test/benchmark_results

cd /home/rsahleanu/rawhash2/test
python3 -u benchmark_segmenters.py

echo "Job finished: $(date)"

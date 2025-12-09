#!/bin/bash
#SBATCH --job-name=blast_validation
#SBATCH --output=blast_validation_%j.log
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general

# BLAST validation workflow for SARS-CoV-2 reads
# Submitted to SLURM for long-running BLAST job

echo "=========================================="
echo "BLAST Validation SLURM Job"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Started: $(date)"
echo ""

# Set input/output
INPUT_SAM="Handley_A2690_54317_Wilen_9_mappings.nonsupplemental.sam"
OUTPUT_PREFIX="Wilen_9_blast"
THREADS=16

# Run the workflow
bash blast_validation_workflow.sh "$INPUT_SAM" "$OUTPUT_PREFIX" "$THREADS"

echo ""
echo "Finished: $(date)"
echo "=========================================="

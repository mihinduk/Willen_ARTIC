#!/bin/bash
#SBATCH --job-name=mmseqs_validation
#SBATCH --output=mmseqs_validation_%j.log
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general

# MMseqs2 validation workflow for SARS-CoV-2 reads
# Usage: sbatch submit_mmseqs_validation.sh input.sam output_prefix [threads]
# MMseqs2 is typically 100-1000x faster than BLAST

# Parse arguments
INPUT_SAM=${1:-"Handley_A2690_54317_Wilen_9_mappings.nonsupplemental.sam"}
OUTPUT_PREFIX=${2:-"Wilen_9_mmseqs"}
THREADS=${3:-16}

echo "=========================================="
echo "MMseqs2 Validation SLURM Job"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Started: $(date)"
echo ""
echo "Input SAM: $INPUT_SAM"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Threads: $THREADS"
echo ""

# Activate MMseqs2 environment
source /ref/sahlab/software/anaconda3/bin/activate
conda activate /ref/sahlab/software/miniforge3/envs/mmseqs2_v15-6f452

# Verify MMseqs2
echo "MMseqs2 version:"
mmseqs version
echo ""

# Run the workflow
bash mmseqs_validation_workflow.sh "$INPUT_SAM" "$OUTPUT_PREFIX" "$THREADS"

echo ""
echo "Finished: $(date)"
echo "=========================================="

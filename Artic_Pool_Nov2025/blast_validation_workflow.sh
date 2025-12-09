#!/bin/bash
#
# Complete BLAST validation workflow for post-ivar trimmed reads
# Usage: bash blast_validation_workflow.sh input.sam output_prefix [threads]
#

if [ $# -lt 2 ]; then
    echo "Usage: bash blast_validation_workflow.sh input.sam output_prefix [threads]"
    echo ""
    echo "Arguments:"
    echo "  input.sam      - Input SAM file (post-ivar trimming)"
    echo "  output_prefix  - Prefix for all output files"
    echo "  threads        - Number of threads for BLAST (default: 8)"
    echo ""
    echo "Example:"
    echo "  bash blast_validation_workflow.sh Handley_A2690_54317_Wilen_9_mappings.final.sam Wilen_9_blast 16"
    echo ""
    echo "This script will:"
    echo "  1. Extract reads from SAM as FASTA"
    echo "  2. Run BLAST against local nt database"
    echo "  3. Parse results and classify hits"
    echo "  4. Generate COVID gene mapping TSV"
    echo "  5. Generate non-COVID category FASTAs"
    exit 1
fi

INPUT_SAM=$1
OUTPUT_PREFIX=$2
THREADS=${3:-8}  # Default to 8 threads if not specified

echo "=========================================="
echo "BLAST Validation Workflow"
echo "=========================================="
echo "Input SAM: $INPUT_SAM"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Threads: $THREADS"
echo ""

# Check if input exists
if [ ! -f "$INPUT_SAM" ]; then
    echo "ERROR: Input SAM file not found: $INPUT_SAM"
    exit 1
fi

# Step 1: Extract reads
echo "=========================================="
echo "Step 1: Extracting reads from SAM..."
echo "=========================================="
echo ""

READS_FASTA="${OUTPUT_PREFIX}_reads.fasta"
METADATA_TSV="${OUTPUT_PREFIX}_metadata.tsv"

python3 extract_reads_for_blast.py "$INPUT_SAM" "$READS_FASTA" \
    --min-mapq 30 \
    --metadata-tsv "$METADATA_TSV"

if [ ! -f "$READS_FASTA" ]; then
    echo "ERROR: Failed to extract reads"
    exit 1
fi

echo ""

# Step 2: BLAST against nt
echo "=========================================="
echo "Step 2: BLASTing against local nt database..."
echo "=========================================="
echo ""

BLAST_RESULTS="${OUTPUT_PREFIX}_blast_results.txt"
NT_DB="/ref/sahlab/data/nt/nt"

# Check if local nt database exists (check for .nal or .00.nhr files)
if [ ! -f "${NT_DB}.nal" ] && [ ! -f "${NT_DB}.00.nhr" ]; then
    echo "WARNING: Local nt database not found at: $NT_DB"
    echo "Falling back to remote BLAST (will be very slow)..."
    blastn -query "$READS_FASTA" \
        -db nt \
        -remote \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \
        -max_target_seqs 5 \
        -out "$BLAST_RESULTS"
else
    echo "Using local nt database: $NT_DB"
    echo "Threads: $THREADS"
    echo ""
    blastn -query "$READS_FASTA" \
        -db "$NT_DB" \
        -num_threads "$THREADS" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \
        -max_target_seqs 5 \
        -out "$BLAST_RESULTS"
fi

if [ ! -f "$BLAST_RESULTS" ]; then
    echo "ERROR: BLAST failed"
    exit 1
fi

echo ""
echo "BLAST complete!"
echo ""

# Step 3: Parse BLAST results
echo "=========================================="
echo "Step 3: Parsing BLAST results..."
echo "=========================================="
echo ""

python3 parse_blast_results.py "$BLAST_RESULTS" "$READS_FASTA" "$METADATA_TSV" "$OUTPUT_PREFIX"

echo ""
echo "=========================================="
echo "WORKFLOW COMPLETE!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  ${OUTPUT_PREFIX}_covid_gene_mapping.tsv    - COVID reads with gene mapping"
echo "  ${OUTPUT_PREFIX}_human_reads.fasta         - Human-matching reads (if any)"
echo "  ${OUTPUT_PREFIX}_bacterial_reads.fasta     - Bacterial-matching reads (if any)"
echo "  ${OUTPUT_PREFIX}_viral_reads.fasta         - Other viral reads (if any)"
echo "  ${OUTPUT_PREFIX}_other_reads.fasta         - Other organism reads (if any)"
echo ""

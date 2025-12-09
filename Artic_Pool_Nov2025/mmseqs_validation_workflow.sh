#!/bin/bash
#
# MMseqs2 validation workflow for post-ivar trimmed reads
# Usage: bash mmseqs_validation_workflow.sh input.sam output_prefix [threads]
#

if [ $# -lt 2 ]; then
    echo "Usage: bash mmseqs_validation_workflow.sh input.sam output_prefix [threads]"
    echo ""
    echo "Arguments:"
    echo "  input.sam      - Input SAM file (post-ivar trimming)"
    echo "  output_prefix  - Prefix for all output files"
    echo "  threads        - Number of threads for MMseqs2 (default: 16)"
    echo ""
    echo "Example:"
    echo "  bash mmseqs_validation_workflow.sh Handley_A2690_54317_Wilen_9_mappings.nonsupplemental.sam Wilen_9_mmseqs 16"
    echo ""
    echo "This script will:"
    echo "  1. Extract reads from SAM as FASTA"
    echo "  2. Create MMseqs2 query database"
    echo "  3. Search against nt database using MMseqs2"
    echo "  4. Parse results and classify hits"
    echo "  5. Generate COVID gene mapping TSV"
    echo "  6. Generate non-COVID category FASTAs"
    exit 1
fi

INPUT_SAM=$1
OUTPUT_PREFIX=$2
THREADS=${3:-16}

echo "=========================================="
echo "MMseqs2 Validation Workflow"
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

# Check if mmseqs is available
if ! command -v mmseqs &> /dev/null; then
    echo "ERROR: mmseqs not found in PATH"
    echo "Please install MMseqs2 or load the appropriate module"
    echo "  conda install -c bioconda mmseqs2"
    echo "  OR module load mmseqs2"
    exit 1
fi

echo "MMseqs2 version:"
mmseqs version
echo ""

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

# Step 2: Create MMseqs2 query database
echo "=========================================="
echo "Step 2: Creating MMseqs2 query database..."
echo "=========================================="
echo ""

QUERY_DB="${OUTPUT_PREFIX}_queryDB"
mmseqs createdb "$READS_FASTA" "$QUERY_DB"

if [ ! -f "${QUERY_DB}.dbtype" ]; then
    echo "ERROR: Failed to create query database"
    exit 1
fi

echo ""

# Step 3: Search against nt database
echo "=========================================="
echo "Step 3: Searching against nt database..."
echo "=========================================="
echo ""

# Check for MMseqs2 nt database
NT_MMSEQS_DB="/ref/sahlab/data/nt/mmseqs_nt_db/nt"
NT_FASTA="/ref/sahlab/data/nt/nt.gz"

if [ -f "${NT_MMSEQS_DB}.dbtype" ]; then
    echo "Using existing MMseqs2 nt database: $NT_MMSEQS_DB"
    TARGET_DB="$NT_MMSEQS_DB"
elif [ -f "$NT_FASTA" ]; then
    echo "MMseqs2 nt database not found."
    echo "Creating MMseqs2 database from: $NT_FASTA"
    echo "This will take 1-2 hours and requires ~500GB disk space..."
    echo ""

    TARGET_DB="${OUTPUT_PREFIX}_nt_mmseqs_db"
    gunzip -c "$NT_FASTA" | mmseqs createdb stdin "$TARGET_DB" --dbtype 2

    if [ ! -f "${TARGET_DB}.dbtype" ]; then
        echo "ERROR: Failed to create target database"
        exit 1
    fi
else
    echo "ERROR: No nt database found!"
    echo "Looked for:"
    echo "  - MMseqs2 database: ${NT_MMSEQS_DB}.dbtype"
    echo "  - FASTA file: $NT_FASTA"
    exit 1
fi

echo ""
echo "Running MMseqs2 search..."
RESULT_DB="${OUTPUT_PREFIX}_resultDB"
TMP_DIR="${OUTPUT_PREFIX}_tmp"

mmseqs search "$QUERY_DB" "$TARGET_DB" "$RESULT_DB" "$TMP_DIR" \
    --threads "$THREADS" \
    --max-seqs 5 \
    -s 7.5 \
    -e 0.001 \
    --search-type 3

if [ ! -f "${RESULT_DB}.dbtype" ]; then
    echo "ERROR: MMseqs2 search failed"
    exit 1
fi

echo ""

# Step 4: Convert to tabular format with taxonomy
echo "=========================================="
echo "Step 4: Converting results to tabular format..."
echo "=========================================="
echo ""

MMSEQS_RESULTS="${OUTPUT_PREFIX}_mmseqs_results.m8"

# Convert to BLAST-like format
mmseqs convertalis "$QUERY_DB" "$TARGET_DB" "$RESULT_DB" "$MMSEQS_RESULTS" \
    --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits" \
    --threads "$THREADS"

if [ ! -f "$MMSEQS_RESULTS" ]; then
    echo "ERROR: Failed to convert results"
    exit 1
fi

echo ""
echo "MMseqs2 search complete!"
echo "Results: $MMSEQS_RESULTS"
echo ""

# Step 5: Add taxonomy information
echo "=========================================="
echo "Step 5: Adding taxonomy information..."
echo "=========================================="
echo ""

MMSEQS_WITH_TAX="${OUTPUT_PREFIX}_mmseqs_with_taxonomy.txt"

python3 add_taxonomy_to_mmseqs.py "$MMSEQS_RESULTS" "$MMSEQS_WITH_TAX"

if [ ! -f "$MMSEQS_WITH_TAX" ]; then
    echo "WARNING: Taxonomy annotation failed, using results without taxonomy"
    MMSEQS_WITH_TAX="$MMSEQS_RESULTS"
fi

echo ""

# Step 6: Parse results
echo "=========================================="
echo "Step 6: Parsing results and classifying hits..."
echo "=========================================="
echo ""

python3 parse_mmseqs_results.py "$MMSEQS_WITH_TAX" "$READS_FASTA" "$METADATA_TSV" "$OUTPUT_PREFIX"

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
echo "MMseqs2 intermediate files:"
echo "  ${OUTPUT_PREFIX}_queryDB.*"
echo "  ${OUTPUT_PREFIX}_resultDB.*"
echo "  ${OUTPUT_PREFIX}_tmp/ (can be deleted)"
echo ""

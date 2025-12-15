#!/bin/bash
#
# Extract reads with high-quality COVID hits into a FASTA file
#

if [ $# -lt 3 ]; then
    echo "Usage: bash extract_covid_reads.sh taxonomy_file.txt reads.fasta output.fasta"
    echo ""
    echo "Example:"
    echo "  bash extract_covid_reads.sh Wilen_6_mmseqs_mmseqs_with_taxonomy_FIXED.txt \\"
    echo "                               Wilen_6_mmseqs_reads.fasta \\"
    echo "                               Wilen_6_covid_hq.fasta"
    exit 1
fi

TAXONOMY_FILE=$1
READS_FASTA=$2
OUTPUT_FASTA=$3

ORGANISM="Severe acute respiratory syndrome coronavirus 2"
MIN_LENGTH=50
MIN_PIDENT=90.0

echo "Extracting high-quality COVID reads..."
echo "  Input taxonomy: $TAXONOMY_FILE"
echo "  Input reads: $READS_FASTA"
echo "  Output: $OUTPUT_FASTA"
echo "  Filters: length >= ${MIN_LENGTH}bp, pident >= ${MIN_PIDENT}%"
echo ""

# Step 1: Extract unique query IDs with high-quality COVID hits
TEMP_IDS=$(mktemp)

grep "$ORGANISM" "$TAXONOMY_FILE" | \
  awk -F'\t' -v len="$MIN_LENGTH" -v pid="$MIN_PIDENT" \
    '$3 >= pid && $4 >= len {print $1}' | \
  sort -u > "$TEMP_IDS"

TOTAL_IDS=$(wc -l < "$TEMP_IDS")
echo "Found $TOTAL_IDS unique reads with high-quality COVID hits"

if [ $TOTAL_IDS -eq 0 ]; then
    echo "No high-quality COVID reads found!"
    rm "$TEMP_IDS"
    exit 1
fi

echo ""
echo "Extracting sequences from FASTA..."

# Step 2: Extract sequences using seqtk (if available) or awk
if command -v seqtk &> /dev/null; then
    seqtk subseq "$READS_FASTA" "$TEMP_IDS" > "$OUTPUT_FASTA"
else
    # Fallback: use awk
    awk 'NR==FNR {ids[$1]=1; next} 
         /^>/ {p=0; id=substr($1,2); if (ids[id]) p=1} 
         p' "$TEMP_IDS" "$READS_FASTA" > "$OUTPUT_FASTA"
fi

# Count sequences in output
SEQS_EXTRACTED=$(grep -c "^>" "$OUTPUT_FASTA")

echo "Extracted $SEQS_EXTRACTED sequences to: $OUTPUT_FASTA"

# Cleanup
rm "$TEMP_IDS"

echo ""
echo "Done!"

#!/bin/bash
#
# Test different minimum mapped length thresholds
# Usage: bash test_length_filters.sh input.sam
#

if [ $# -ne 1 ]; then
    echo "Usage: bash test_length_filters.sh input.sam"
    echo ""
    echo "Example:"
    echo "  bash test_length_filters.sh Handley_A2690_54317_Wilen_9_mappings.nonsupplemental.sam"
    exit 1
fi

INPUT_SAM=$1
BASENAME=$(basename "$INPUT_SAM" .sam)

echo "=========================================="
echo "Testing Length Filter Thresholds"
echo "=========================================="
echo "Input: $INPUT_SAM"
echo ""

# Test thresholds
THRESHOLDS=(150 175 200 250)

# Step 1: Filter SAM files
echo "Step 1: Filtering SAM files..."
echo ""

for MIN_LEN in "${THRESHOLDS[@]}"; do
    OUTPUT_SAM="${BASENAME}_min${MIN_LEN}.sam"

    echo "  Filtering ≥${MIN_LEN}bp → ${OUTPUT_SAM}"
    python3 filter_sam_by_length.py "$INPUT_SAM" "$MIN_LEN" "$OUTPUT_SAM"
    echo ""
done

# Step 2: Generate coverage maps
echo ""
echo "=========================================="
echo "Step 2: Generating coverage maps..."
echo "=========================================="
echo ""

for MIN_LEN in "${THRESHOLDS[@]}"; do
    FILTERED_SAM="${BASENAME}_min${MIN_LEN}.sam"

    if [ -f "$FILTERED_SAM" ]; then
        echo "  Plotting coverage for min ${MIN_LEN}bp..."
        python3 plot_coverage_map.py "$FILTERED_SAM" \
            --accession MN908947.3 \
            --min-mapq 30
        echo ""
    else
        echo "  ERROR: $FILTERED_SAM not found!"
    fi
done

# Step 3: Calculate median mapped lengths
echo ""
echo "=========================================="
echo "Step 3: Calculating median mapped lengths"
echo "=========================================="
echo ""

for MIN_LEN in "${THRESHOLDS[@]}"; do
    FILTERED_SAM="${BASENAME}_min${MIN_LEN}.sam"

    if [ -f "$FILTERED_SAM" ]; then
        echo "  Min ${MIN_LEN}bp threshold:"
        python3 calculate_median_mapped_length.py < "$FILTERED_SAM" | grep -E "Reads analyzed|Median|Mean"
        echo ""
    fi
done

# Summary
echo ""
echo "=========================================="
echo "COMPLETE!"
echo "=========================================="
echo ""
echo "Generated files:"
for MIN_LEN in "${THRESHOLDS[@]}"; do
    echo "  ${BASENAME}_min${MIN_LEN}.sam"
    echo "  ${BASENAME}_min${MIN_LEN}_coverage.png"
    echo "  ${BASENAME}_min${MIN_LEN}_coverage.pdf"
done
echo ""
echo "Compare coverage maps to see effect of filtering!"

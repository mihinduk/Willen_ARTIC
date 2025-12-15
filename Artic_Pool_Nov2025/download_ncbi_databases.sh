#!/bin/bash
#
# Download NCBI taxonomy databases for annotation
# This only needs to be run once (databases can be shared across projects)
#

TAXONOMY_DIR="/scratch/sahlab/kathie/ncbi_taxonomy"
mkdir -p "$TAXONOMY_DIR"
cd "$TAXONOMY_DIR"

echo "=========================================="
echo "Downloading NCBI Taxonomy Databases"
echo "=========================================="
echo "Target directory: $TAXONOMY_DIR"
echo ""

# Download accession2taxid (nucleotide GenBank)
echo "Step 1: Downloading nucl_gb.accession2taxid.gz (~11 GB)..."
echo "  This contains mapping of GenBank accessions to taxonomy IDs"
if [ ! -f nucl_gb.accession2taxid.gz ]; then
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    echo "  Downloaded successfully"
else
    echo "  File already exists, skipping download"
fi

# Download taxonomy dump
echo ""
echo "Step 2: Downloading NCBI taxonomy dump (~50 MB)..."
echo "  This contains taxid to organism name mappings"
if [ ! -f taxdmp.zip ]; then
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
    echo "  Downloaded successfully"
else
    echo "  File already exists, skipping download"
fi

# Extract names.dmp
echo ""
echo "Step 3: Extracting names.dmp..."
if [ ! -f names.dmp ]; then
    unzip -o taxdmp.zip names.dmp
    echo "  Extracted successfully"
else
    echo "  File already exists, skipping extraction"
fi

echo ""
echo "=========================================="
echo "DOWNLOAD COMPLETE"
echo "=========================================="
echo ""
echo "Files created:"
ls -lh nucl_gb.accession2taxid.gz names.dmp taxdmp.zip 2>/dev/null
echo ""
echo "You can now use add_ncbi_taxonomy.py with these databases:"
echo ""
echo "python3 add_ncbi_taxonomy.py results.m8 output_with_tax.txt \\"
echo "    --accession2taxid $TAXONOMY_DIR/nucl_gb.accession2taxid.gz \\"
echo "    --names $TAXONOMY_DIR/names.dmp"
echo ""

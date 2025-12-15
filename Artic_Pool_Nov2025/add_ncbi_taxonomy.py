#!/usr/bin/env python3
"""
Add proper NCBI taxonomy to MMseqs2 results using accession2taxid database.

This script:
1. Parses MMseqs2 m8 results
2. Looks up accessions in NCBI's nucl_gb.accession2taxid database
3. Gets organism names from NCBI taxonomy names.dmp
4. Outputs BLAST-compatible format with taxonomy

Usage:
    python3 add_ncbi_taxonomy.py results.m8 output_with_tax.txt \\
        --accession2taxid /path/to/nucl_gb.accession2taxid.gz \\
        --names /path/to/names.dmp

Dependencies:
    Download these files:
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
    unzip taxdmp.zip names.dmp
"""

import sys
import gzip
import argparse
from collections import defaultdict

def load_names_dmp(names_file):
    """Load taxid -> scientific name mapping from names.dmp."""
    print(f"Loading taxonomy names from {names_file}...")
    taxid_to_name = {}
    
    with open(names_file, 'r') as f:
        for line in f:
            fields = [x.strip() for x in line.split('|')]
            if len(fields) >= 4 and fields[3] == 'scientific name':
                taxid = fields[0]
                name = fields[1]
                taxid_to_name[taxid] = name
    
    print(f"  Loaded {len(taxid_to_name):,} taxonomy names")
    return taxid_to_name

def build_accession_index(accession2taxid_file, needed_accessions):
    """
    Build index of accession -> taxid for only the accessions we need.
    This is much faster than loading the entire 40GB+ file.
    """
    print(f"\nBuilding accession index from {accession2taxid_file}...")
    print(f"  Looking for {len(needed_accessions):,} unique accessions...")
    
    accession_to_taxid = {}
    found = 0
    total_lines = 0
    
    # Open with gzip if it's compressed
    if accession2taxid_file.endswith('.gz'):
        f = gzip.open(accession2taxid_file, 'rt')
    else:
        f = open(accession2taxid_file, 'r')
    
    # Skip header
    next(f)
    
    for line in f:
        total_lines += 1
        fields = line.strip().split('\t')
        if len(fields) < 3:
            continue
        
        accession = fields[0]
        accession_version = fields[1]  # With version number
        taxid = fields[2]
        
        # Check if this is one of our needed accessions
        # Try both with and without version number
        if accession in needed_accessions or accession_version in needed_accessions:
            accession_to_taxid[accession] = taxid
            accession_to_taxid[accession_version] = taxid
            found += 1
        
        # Progress update
        if total_lines % 10000000 == 0:
            print(f"    Scanned {total_lines:,} lines, found {found:,}/{len(needed_accessions):,} accessions...")
        
        # Stop early if we found everything
        if found >= len(needed_accessions):
            print(f"    Found all needed accessions!")
            break
    
    f.close()
    
    print(f"  Scanned {total_lines:,} total lines")
    print(f"  Found {found:,}/{len(needed_accessions):,} accessions ({100*found/len(needed_accessions):.1f}%)")
    
    return accession_to_taxid

def parse_accession(target_id):
    """Extract clean accession from target ID."""
    # Remove any prefixes and extract the accession
    # Format examples: OY655439.1, gi|12345|ref|NC_045512.2|
    import re
    match = re.search(r'([A-Z]{1,2}_?\d+\.?\d*)', target_id)
    if match:
        return match.group(1)
    return target_id.split('|')[0].split()[0]

def main():
    parser = argparse.ArgumentParser(description='Add NCBI taxonomy to MMseqs2 results')
    parser.add_argument('input_m8', help='Input MMseqs2 m8 format file')
    parser.add_argument('output', help='Output file with taxonomy')
    parser.add_argument('--accession2taxid', required=True, 
                       help='Path to nucl_gb.accession2taxid.gz')
    parser.add_argument('--names', required=True,
                       help='Path to names.dmp from NCBI taxonomy')
    
    args = parser.parse_args()
    
    print("="*70)
    print("NCBI Taxonomy Annotation for MMseqs2 Results")
    print("="*70)
    
    # Step 1: Load taxonomy names
    taxid_to_name = load_names_dmp(args.names)
    
    # Step 2: First pass - collect all unique accessions we need to look up
    print(f"\nFirst pass: collecting accessions from {args.input_m8}...")
    needed_accessions = set()
    total_hits = 0
    
    with open(args.input_m8, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                total_hits += 1
                target = fields[1]
                accession = parse_accession(target)
                needed_accessions.add(accession)
    
    print(f"  Found {total_hits:,} total hits")
    print(f"  Found {len(needed_accessions):,} unique accessions to look up")
    
    # Step 3: Build accession index (only for accessions we need)
    accession_to_taxid = build_accession_index(args.accession2taxid, needed_accessions)
    
    # Step 4: Second pass - annotate results
    print(f"\nSecond pass: annotating results...")
    annotated = 0
    unknown = 0
    
    organism_counts = defaultdict(int)
    
    with open(args.input_m8, 'r') as infile, open(args.output, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            
            # Parse MMseqs2 fields
            query, target = fields[0], fields[1]
            pident, alnlen, mismatch, gapopen = fields[2:6]
            qstart, qend, tstart, tend = fields[6:10]
            evalue, bits = fields[10:12]
            
            # Get accession and look up taxonomy
            accession = parse_accession(target)
            taxid = accession_to_taxid.get(accession, '0')
            organism = taxid_to_name.get(taxid, 'Unknown organism')
            
            # Track stats
            organism_counts[organism] += 1
            if taxid != '0':
                annotated += 1
            else:
                unknown += 1
            
            # Write BLAST-compatible format
            outfile.write(f"{query}\t{target}\t{pident}\t{alnlen}\t{mismatch}\t{gapopen}\t"
                         f"{qstart}\t{qend}\t{tstart}\t{tend}\t{evalue}\t{bits}\t"
                         f"{taxid}\t{organism}\n")
    
    print(f"\n{'='*70}")
    print(f"ANNOTATION COMPLETE")
    print(f"{'='*70}")
    print(f"Total hits: {total_hits:,}")
    print(f"Successfully annotated: {annotated:,} ({100*annotated/total_hits:.1f}%)")
    print(f"Unknown taxonomy: {unknown:,} ({100*unknown/total_hits:.1f}%)")
    
    print(f"\nTop 10 organisms:")
    for organism, count in sorted(organism_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
        pct = 100 * count / total_hits
        print(f"  {organism}: {count:,} ({pct:.1f}%)")
    
    print(f"\nOutput written to: {args.output}")

if __name__ == '__main__':
    main()

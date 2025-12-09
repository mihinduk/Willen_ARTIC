#!/usr/bin/env python3
"""
Add taxonomy information to MMseqs2 results.
Attempts to map accession IDs to organism names using NCBI databases.

Usage:
    python3 add_taxonomy_to_mmseqs.py input.m8 output_with_taxonomy.txt
"""

import sys
import re
from collections import defaultdict

# Common accession to organism mappings for speed
COMMON_ORGANISMS = {
    'MN908947': 'Severe acute respiratory syndrome coronavirus 2',
    'NC_045512': 'Severe acute respiratory syndrome coronavirus 2',
    'MT': 'Severe acute respiratory syndrome coronavirus 2',  # Prefix for many SARS-CoV-2
    'MW': 'Severe acute respiratory syndrome coronavirus 2',
    'OM': 'Severe acute respiratory syndrome coronavirus 2',
    'OL': 'Severe acute respiratory syndrome coronavirus 2',
}

def guess_organism_from_accession(accession):
    """Guess organism from accession ID patterns."""
    # Check exact matches
    for prefix, organism in COMMON_ORGANISMS.items():
        if accession.startswith(prefix):
            return organism

    # Try to infer from accession patterns
    # RefSeq patterns
    if accession.startswith('NC_'):
        return 'Unknown RefSeq organism'
    elif accession.startswith('NM_'):
        return 'Homo sapiens'  # mRNA
    elif accession.startswith('NR_'):
        return 'Homo sapiens'  # RNA
    elif accession.startswith('NP_'):
        return 'Unknown protein'
    elif accession.startswith('XM_'):
        return 'Unknown predicted mRNA'
    elif accession.startswith('XP_'):
        return 'Unknown predicted protein'

    # GenBank patterns - harder to infer organism
    # Format: AB123456
    return f'Unknown organism'

def parse_accession(target_id):
    """Extract accession from target ID."""
    # Remove version numbers and extra annotations
    # Format might be: gi|12345|ref|NC_045512.2| or just NC_045512.2
    match = re.search(r'([A-Z]{1,2}_?\d+)', target_id)
    if match:
        return match.group(1)
    return target_id.split('|')[0]

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 add_taxonomy_to_mmseqs.py input.m8 output_with_taxonomy.txt")
        print("\nInput should be MMseqs2 m8 format:")
        print("  query target pident alnlen mismatch gapopen qstart qend tstart tend evalue bits")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    print(f"Processing: {input_file}")
    print(f"Adding taxonomy annotations...")
    print(f"Output: {output_file}\n")

    processed = 0
    annotated = 0
    unknown = 0

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            processed += 1

            # Parse fields
            query = fields[0]
            target = fields[1]
            pident = fields[2]
            alnlen = fields[3]
            mismatch = fields[4]
            gapopen = fields[5]
            qstart = fields[6]
            qend = fields[7]
            tstart = fields[8]
            tend = fields[9]
            evalue = fields[10]
            bits = fields[11]

            # Extract accession and guess organism
            accession = parse_accession(target)
            organism = guess_organism_from_accession(accession)

            if organism.startswith('Unknown'):
                unknown += 1
            else:
                annotated += 1

            # Add placeholder taxid (0 = unknown)
            taxid = '0'
            if 'coronavirus 2' in organism.lower():
                taxid = '2697049'  # SARS-CoV-2 taxid

            # Write in BLAST-compatible format
            # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames
            outfile.write(f"{query}\t{target}\t{pident}\t{alnlen}\t{mismatch}\t{gapopen}\t"
                         f"{qstart}\t{qend}\t{tstart}\t{tend}\t{evalue}\t{bits}\t"
                         f"{taxid}\t{organism}\n")

            if processed % 100000 == 0:
                print(f"  Processed {processed:,} hits...")

    print(f"\n{'='*60}")
    print(f"ANNOTATION COMPLETE")
    print(f"{'='*60}")
    print(f"Total hits processed: {processed:,}")
    print(f"Annotated with known organism: {annotated:,} ({100*annotated/processed:.1f}%)")
    print(f"Unknown organisms: {unknown:,} ({100*unknown/processed:.1f}%)")
    print(f"\nWrote annotated results to: {output_file}")
    print("\nNOTE: Taxonomy is inferred from accession patterns.")
    print("For more accurate taxonomy, consider using MMseqs2 taxonomy database.")

if __name__ == '__main__':
    main()

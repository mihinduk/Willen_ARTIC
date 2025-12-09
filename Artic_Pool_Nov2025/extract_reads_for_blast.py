#!/usr/bin/env python3
"""
Extract mapped reads from SAM file for BLAST validation.
Preserves mapping coordinates for downstream gene mapping.

Usage:
    python3 extract_reads_for_blast.py input.sam output.fasta [--min-mapq 30]
"""

import sys
import re
import argparse

def parse_cigar_length(cigar):
    """Calculate mapped length from CIGAR."""
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    return sum(int(num) for num, op in ops if op in ['M', '=', 'X'])

def parse_cigar_end_pos(cigar, start_pos):
    """Calculate end position on reference from CIGAR."""
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    current_pos = start_pos
    for length, op in ops:
        length = int(length)
        if op in ['M', '=', 'X', 'D', 'N']:  # Consumes reference
            current_pos += length
    return current_pos - 1  # Last covered position

def main():
    parser = argparse.ArgumentParser(
        description='Extract mapped reads from SAM for BLAST validation'
    )
    parser.add_argument('sam_file', help='Input SAM file')
    parser.add_argument('output_fasta', help='Output FASTA file')
    parser.add_argument('--min-mapq', type=int, default=30,
                        help='Minimum MAPQ for filtering reads (default: 30)')
    parser.add_argument('--metadata-tsv', help='Optional: output TSV with read metadata')

    args = parser.parse_args()

    print(f"Processing: {args.sam_file}")
    print(f"Min MAPQ: {args.min_mapq}\n")

    total_reads = 0
    passed_reads = 0

    fasta_out = open(args.output_fasta, 'w')
    metadata_out = open(args.metadata_tsv, 'w') if args.metadata_tsv else None

    if metadata_out:
        metadata_out.write("read_id\tstart\tend\tmapped_length\tmapq\n")

    with open(args.sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            total_reads += 1

            read_name = fields[0]
            flag = int(fields[1])
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]

            # Skip unmapped or low quality
            if flag & 0x4 or mapq < args.min_mapq:
                continue

            passed_reads += 1

            # Calculate mapping info
            mapped_length = parse_cigar_length(cigar)
            end_pos = parse_cigar_end_pos(cigar, pos)

            # Create unique read ID with mapping info
            read_id = f"{read_name}_pos{pos}-{end_pos}_len{mapped_length}_mapq{mapq}"

            # Write FASTA
            fasta_out.write(f">{read_id}\n")
            fasta_out.write(f"{seq}\n")

            # Write metadata
            if metadata_out:
                metadata_out.write(f"{read_name}\t{pos}\t{end_pos}\t{mapped_length}\t{mapq}\n")

            if passed_reads % 100000 == 0:
                print(f"  Processed {total_reads:,} reads, wrote {passed_reads:,}...")

    fasta_out.close()
    if metadata_out:
        metadata_out.close()

    print(f"\n{'='*60}")
    print(f"EXTRACTION COMPLETE")
    print(f"{'='*60}")
    print(f"Total reads: {total_reads:,}")
    if total_reads > 0:
        print(f"Passed filters (MAPQ ≥{args.min_mapq}): {passed_reads:,} ({100*passed_reads/total_reads:.2f}%)")
    else:
        print(f"Passed filters (MAPQ ≥{args.min_mapq}): {passed_reads:,}")
        print("\nWARNING: No reads found in SAM file!")
        print("Check that the file contains alignment records (not just headers)")
    print(f"\nWrote {passed_reads:,} sequences to: {args.output_fasta}")
    if args.metadata_tsv:
        print(f"Wrote metadata to: {args.metadata_tsv}")
    if passed_reads > 0:
        print("\nReady for BLAST against nt database!")
    else:
        print("\nERROR: No reads to BLAST!")

if __name__ == '__main__':
    main()

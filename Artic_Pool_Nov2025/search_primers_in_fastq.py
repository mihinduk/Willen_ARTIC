#!/usr/bin/env python3
"""
Search for primer sequences in FASTQ files and report counts.

Usage:
    python3 search_primers_in_fastq.py primers.fasta reads.fastq.gz

Allows up to 1 mismatch per primer.
"""

import sys
import gzip
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def hamming_distance(s1, s2):
    """Calculate Hamming distance between two strings."""
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def search_primer(read_seq, primer_seq, max_mismatches=1):
    """Search for primer in read allowing mismatches."""
    primer_len = len(primer_seq)
    read_len = len(read_seq)

    # Check both orientations
    for seq in [primer_seq, str(Seq(primer_seq).reverse_complement())]:
        # Slide window across read
        for i in range(read_len - primer_len + 1):
            window = read_seq[i:i+primer_len]
            if hamming_distance(window, seq) <= max_mismatches:
                return True

    return False

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 search_primers_in_fastq.py primers.fasta reads.fastq[.gz]", file=sys.stderr)
        sys.exit(1)

    primer_file = sys.argv[1]
    fastq_file = sys.argv[2]

    # Load primers
    print(f"Loading primers from {primer_file}...", file=sys.stderr)
    primers = {}
    for record in SeqIO.parse(primer_file, 'fasta'):
        primers[record.id] = str(record.seq).upper()

    print(f"Loaded {len(primers)} primers", file=sys.stderr)

    # Count primer matches
    primer_counts = defaultdict(int)
    total_reads = 0
    reads_with_primers = 0

    print(f"Searching in {fastq_file}...", file=sys.stderr)

    # Handle gzipped or plain FASTQ
    if fastq_file.endswith('.gz'):
        handle = gzip.open(fastq_file, 'rt')
    else:
        handle = open(fastq_file, 'r')

    try:
        for record in SeqIO.parse(handle, 'fastq'):
            total_reads += 1
            read_seq = str(record.seq).upper()

            found_any = False
            for primer_name, primer_seq in primers.items():
                if search_primer(read_seq, primer_seq, max_mismatches=1):
                    primer_counts[primer_name] += 1
                    found_any = True

            if found_any:
                reads_with_primers += 1

            # Progress report
            if total_reads % 100000 == 0:
                print(f"  Processed {total_reads:,} reads...", file=sys.stderr)

    finally:
        handle.close()

    # Report results
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"RESULTS", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Total reads: {total_reads:,}", file=sys.stderr)
    print(f"Reads with primers: {reads_with_primers:,} ({100*reads_with_primers/total_reads:.1f}%)", file=sys.stderr)
    print(f"\nPrimer counts (allowing 1 mismatch):", file=sys.stderr)

    # Print header
    print("\nPrimer\tCount\tPercent")

    # Sort by count descending
    for primer_name, count in sorted(primer_counts.items(), key=lambda x: x[1], reverse=True):
        pct = 100 * count / total_reads if total_reads > 0 else 0
        print(f"{primer_name}\t{count}\t{pct:.2f}%")

    # Find missing primers
    missing = set(primers.keys()) - set(primer_counts.keys())
    if missing:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"WARNING: {len(missing)} primers NOT found:", file=sys.stderr)
        for primer in sorted(missing):
            print(f"  - {primer}", file=sys.stderr)

if __name__ == '__main__':
    main()

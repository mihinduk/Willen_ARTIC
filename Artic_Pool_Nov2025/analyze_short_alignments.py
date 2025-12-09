#!/usr/bin/env python3
"""
Analyze short alignments to understand what's being mapped.
Extracts reads with mapped length below threshold.

Usage:
    python3 analyze_short_alignments.py input.sam min_length max_length
"""

import sys
import re
from collections import Counter, defaultdict

def parse_cigar(cigar_string):
    """Parse CIGAR to get mapped length."""
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    mapped_length = sum(int(length) for length, op in ops if op in ['M', '=', 'X'])
    return mapped_length, ops

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 analyze_short_alignments.py input.sam min_length max_length")
        print("Example: python3 analyze_short_alignments.py sample.sam 19 30")
        sys.exit(1)

    sam_file = sys.argv[1]
    min_length = int(sys.argv[2])
    max_length = int(sys.argv[3])

    print(f"Analyzing alignments with mapped length: {min_length}-{max_length}bp")
    print(f"Reading: {sam_file}\n")

    # Statistics
    total_reads = 0
    target_reads = 0
    sequences = []
    mapq_dist = Counter()
    position_dist = defaultdict(int)
    cigar_patterns = Counter()

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            total_reads += 1

            flag = int(fields[1])
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]

            # Skip unmapped
            if flag & 0x4:
                continue

            mapped_length, cigar_ops = parse_cigar(cigar)

            if min_length <= mapped_length <= max_length:
                target_reads += 1
                sequences.append(seq)
                mapq_dist[mapq] += 1
                # Bin positions by 1000bp
                position_bin = (pos // 1000) * 1000
                position_dist[position_bin] += 1
                cigar_patterns[cigar] += 1

                # Show first 10 examples
                if target_reads <= 10:
                    print(f"Example {target_reads}:")
                    print(f"  Position: {pos}")
                    print(f"  MAPQ: {mapq}")
                    print(f"  CIGAR: {cigar}")
                    print(f"  Mapped length: {mapped_length}bp")
                    print(f"  Sequence: {seq[:50]}...")
                    print()

    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Total reads processed: {total_reads:,}")
    print(f"Reads with {min_length}-{max_length}bp mapped: {target_reads:,} ({100*target_reads/total_reads:.2f}%)")

    print(f"\n{'='*60}")
    print(f"MAPQ Distribution")
    print(f"{'='*60}")
    for mapq, count in sorted(mapq_dist.items(), reverse=True):
        print(f"MAPQ {mapq:2d}: {count:6,} reads ({100*count/target_reads:.1f}%)")

    print(f"\n{'='*60}")
    print(f"Genome Position Distribution (1kb bins)")
    print(f"{'='*60}")
    for pos, count in sorted(position_dist.items())[:20]:
        print(f"{pos:5,}-{pos+1000:5,}bp: {count:6,} reads")

    print(f"\n{'='*60}")
    print(f"Top 10 CIGAR Patterns")
    print(f"{'='*60}")
    for cigar, count in cigar_patterns.most_common(10):
        print(f"{cigar:20s}: {count:6,} reads")

    # Sequence analysis
    if sequences:
        print(f"\n{'='*60}")
        print(f"Sequence Content Analysis")
        print(f"{'='*60}")

        # Check for adapter sequences and homopolymers
        adapter_patterns = {
            'Illumina_adapter': 'AGATCGGAAGAGC',
            'Nextera_adapter': 'CTGTCTCTTATACACATCT',
            'Poly-A': 'AAAAAAAAAA',
            'Poly-T': 'TTTTTTTTTT',
            'Poly-G': 'GGGGGGGGGG',
            'Poly-C': 'CCCCCCCCCC'
        }

        poly_g_count = 0
        poly_c_count = 0
        poly_g_or_c_count = 0

        for name, pattern in adapter_patterns.items():
            count = sum(1 for seq in sequences if pattern in seq)
            if name == 'Poly-G':
                poly_g_count = count
            elif name == 'Poly-C':
                poly_c_count = count
            if count > 0:
                print(f"{name}: {count} sequences ({100*count/len(sequences):.1f}%)")

        # Count sequences with poly-G OR poly-C
        poly_g_or_c_count = sum(1 for seq in sequences
                                if 'GGGGGGGGGG' in seq or 'CCCCCCCCCC' in seq)
        print(f"\nPoly-G or Poly-C: {poly_g_or_c_count} sequences ({100*poly_g_or_c_count/len(sequences):.1f}%)")

        # GC content
        total_bases = sum(len(seq) for seq in sequences)
        gc_count = sum(seq.count('G') + seq.count('C') for seq in sequences)
        print(f"GC content: {100*gc_count/total_bases:.1f}%")

if __name__ == '__main__':
    main()

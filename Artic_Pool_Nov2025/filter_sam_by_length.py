#!/usr/bin/env python3
"""
Filter SAM file by minimum mapped length.
Tests effect of different length thresholds on coverage.

Usage:
    python3 filter_sam_by_length.py input.sam min_mapped_length output.sam
"""

import sys
import re

def parse_cigar_length(cigar):
    """Calculate mapped length from CIGAR."""
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    return sum(int(num) for num, op in ops if op in ['M', '=', 'X'])

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 filter_sam_by_length.py input.sam min_mapped_length output.sam")
        print("\nExample: Test 150bp threshold")
        print("  python3 filter_sam_by_length.py sample.sam 150 sample_min150.sam")
        sys.exit(1)

    input_sam = sys.argv[1]
    min_length = int(sys.argv[2])
    output_sam = sys.argv[3]

    print(f"Filtering: {input_sam}")
    print(f"Minimum mapped length: {min_length}bp")
    print(f"Output: {output_sam}\n")

    total_reads = 0
    passed_reads = 0
    failed_unmapped = 0
    failed_length = 0

    with open(input_sam, 'r') as infile, open(output_sam, 'w') as outfile:
        for line in infile:
            # Keep header lines
            if line.startswith('@'):
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            total_reads += 1

            flag = int(fields[1])
            cigar = fields[5]

            # Skip unmapped
            if flag & 0x4:
                failed_unmapped += 1
                continue

            # Check mapped length
            mapped_length = parse_cigar_length(cigar)

            if mapped_length >= min_length:
                outfile.write(line)
                passed_reads += 1
            else:
                failed_length += 1

            if total_reads % 100000 == 0:
                print(f"  Processed {total_reads:,} reads...")

    print(f"\n{'='*60}")
    print(f"FILTERING RESULTS")
    print(f"{'='*60}")
    print(f"Total reads: {total_reads:,}")
    print(f"Unmapped: {failed_unmapped:,} ({100*failed_unmapped/total_reads:.1f}%)")
    print(f"Mapped <{min_length}bp: {failed_length:,} ({100*failed_length/total_reads:.1f}%)")
    print(f"Passed (≥{min_length}bp): {passed_reads:,} ({100*passed_reads/total_reads:.1f}%)")
    print(f"\nWrote {passed_reads:,} reads to: {output_sam}")

if __name__ == '__main__':
    main()

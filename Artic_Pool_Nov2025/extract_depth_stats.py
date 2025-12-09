#!/usr/bin/env python3
"""
Extract depth statistics from SAM files.
Outputs position-by-position depth and summary statistics.

Usage:
    python3 extract_depth_stats.py input.sam --genome-length 29903 [--min-depth 50] [--min-mapq 30]
"""

import sys
import re
import argparse
from collections import defaultdict

def parse_cigar(cigar_string, start_pos):
    """Parse CIGAR string and return list of covered positions."""
    positions = []
    current_pos = start_pos

    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    for length, op in ops:
        length = int(length)
        if op in ['M', '=', 'X']:  # Match/mismatch - covers reference
            positions.extend(range(current_pos, current_pos + length))
            current_pos += length
        elif op in ['D', 'N']:  # Deletion/skip - advances reference but no coverage
            current_pos += length
        # I, S, H, P don't advance reference position

    return positions

def parse_sam_coverage(sam_file, min_mapq=30):
    """Parse SAM file and calculate coverage at each position."""
    coverage = defaultdict(int)
    total_reads = 0
    mapped_reads = 0

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

            # Skip unmapped or low quality
            if flag & 0x4 or mapq < min_mapq:
                continue

            mapped_reads += 1

            # Get all positions covered by this read
            positions = parse_cigar(cigar, pos)
            for p in positions:
                coverage[p] += 1

    return coverage, total_reads, mapped_reads

def calculate_coverage_stats(coverage, genome_length, min_depth):
    """Calculate summary statistics for coverage."""
    # Create full genome array
    depth_array = [coverage.get(i, 0) for i in range(1, genome_length + 1)]

    # Basic stats
    total_bases = genome_length
    covered_bases = sum(1 for d in depth_array if d > 0)
    bases_above_threshold = sum(1 for d in depth_array if d >= min_depth)

    # Mean and median depth
    mean_depth = sum(depth_array) / total_bases if total_bases > 0 else 0
    sorted_depths = sorted(depth_array)
    median_depth = sorted_depths[total_bases // 2] if total_bases > 0 else 0

    # Coverage distribution
    coverage_dist = {}
    for threshold in [1, 10, 20, 30, 50, 100, 200]:
        count = sum(1 for d in depth_array if d >= threshold)
        coverage_dist[threshold] = count

    # Find coverage gaps (regions with 0 coverage)
    gaps = []
    in_gap = False
    gap_start = None

    for pos in range(1, genome_length + 1):
        depth = coverage.get(pos, 0)
        if depth == 0:
            if not in_gap:
                gap_start = pos
                in_gap = True
        else:
            if in_gap:
                gaps.append((gap_start, pos - 1))
                in_gap = False

    # Close final gap if at end of genome
    if in_gap:
        gaps.append((gap_start, genome_length))

    return {
        'total_bases': total_bases,
        'covered_bases': covered_bases,
        'bases_above_threshold': bases_above_threshold,
        'mean_depth': mean_depth,
        'median_depth': median_depth,
        'max_depth': max(depth_array) if depth_array else 0,
        'coverage_dist': coverage_dist,
        'gaps': gaps
    }

def main():
    parser = argparse.ArgumentParser(description='Extract depth statistics from SAM file')
    parser.add_argument('sam_file', help='Input SAM file')
    parser.add_argument('--genome-length', type=int, required=True,
                        help='Length of reference genome (e.g., 29903 for SARS-CoV-2)')
    parser.add_argument('--min-depth', type=int, default=50,
                        help='Minimum depth threshold for coverage stats (default: 50)')
    parser.add_argument('--min-mapq', type=int, default=30,
                        help='Minimum MAPQ for filtering reads (default: 30)')
    parser.add_argument('--output-depths', help='Output file for position-by-position depths (TSV)')
    parser.add_argument('--no-summary', action='store_true',
                        help='Skip summary statistics output')

    args = parser.parse_args()

    print(f"Processing: {args.sam_file}")
    print(f"Genome length: {args.genome_length:,} bp")
    print(f"Min MAPQ: {args.min_mapq}")
    print(f"Coverage threshold: {args.min_depth}X\n")

    # Parse SAM file
    print("Parsing SAM file...")
    coverage, total_reads, mapped_reads = parse_sam_coverage(args.sam_file, args.min_mapq)
    print(f"Total reads: {total_reads:,}")
    print(f"Mapped reads (MAPQ ≥{args.min_mapq}): {mapped_reads:,}\n")

    # Calculate statistics
    stats = calculate_coverage_stats(coverage, args.genome_length, args.min_depth)

    # Output position-by-position depths
    if args.output_depths:
        print(f"Writing position-by-position depths to: {args.output_depths}")
        with open(args.output_depths, 'w') as f:
            f.write("Position\tDepth\n")
            for pos in range(1, args.genome_length + 1):
                depth = coverage.get(pos, 0)
                f.write(f"{pos}\t{depth}\n")
        print(f"Wrote {args.genome_length:,} positions\n")

    # Output summary statistics
    if not args.no_summary:
        print("=" * 70)
        print("COVERAGE SUMMARY")
        print("=" * 70)
        print(f"Total bases in genome: {stats['total_bases']:,}")
        print(f"Bases with coverage > 0: {stats['covered_bases']:,} "
              f"({100 * stats['covered_bases'] / stats['total_bases']:.2f}%)")
        print(f"Bases with coverage ≥ {args.min_depth}X: {stats['bases_above_threshold']:,} "
              f"({100 * stats['bases_above_threshold'] / stats['total_bases']:.2f}%)")
        print()
        print(f"Mean depth: {stats['mean_depth']:.2f}X")
        print(f"Median depth: {stats['median_depth']:.0f}X")
        print(f"Max depth: {stats['max_depth']:,}X")
        print()

        print("=" * 70)
        print("COVERAGE DISTRIBUTION")
        print("=" * 70)
        for threshold in sorted(stats['coverage_dist'].keys()):
            count = stats['coverage_dist'][threshold]
            pct = 100 * count / stats['total_bases']
            print(f"≥{threshold:3d}X: {count:7,} bases ({pct:6.2f}%)")
        print()

        # Coverage gaps
        if stats['gaps']:
            print("=" * 70)
            print(f"COVERAGE GAPS (0X coverage)")
            print("=" * 70)
            print(f"Total gaps: {len(stats['gaps'])}")
            print(f"\nLargest gaps (≥100bp):")
            large_gaps = [(start, end, end - start + 1) for start, end in stats['gaps']
                         if end - start + 1 >= 100]
            large_gaps.sort(key=lambda x: x[2], reverse=True)

            if large_gaps:
                for start, end, length in large_gaps[:20]:  # Show top 20
                    print(f"  {start:7,} - {end:7,} ({length:6,} bp)")
            else:
                print("  None ≥100bp")
            print()

            # Summary of gap sizes
            gap_sizes = [end - start + 1 for start, end in stats['gaps']]
            total_gap_bases = sum(gap_sizes)
            print(f"Total bases in gaps: {total_gap_bases:,} "
                  f"({100 * total_gap_bases / stats['total_bases']:.2f}%)")
        else:
            print("=" * 70)
            print("No coverage gaps found!")
            print("=" * 70)

if __name__ == '__main__':
    main()

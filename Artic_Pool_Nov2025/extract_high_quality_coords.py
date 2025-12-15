#!/usr/bin/env python3
"""
Extract and analyze coordinates of high-quality hits for a specific organism.

Usage:
    python3 extract_high_quality_coords.py taxonomy_file.txt "organism_name" [output.tsv]
"""

import sys
from collections import defaultdict

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 extract_high_quality_coords.py taxonomy_file.txt \"organism_name\" [output.tsv]")
        print("\nExample:")
        print("  python3 extract_high_quality_coords.py results.txt \"Severe acute respiratory syndrome coronavirus 2\"")
        sys.exit(1)
    
    input_file = sys.argv[1]
    organism = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    # Thresholds
    MIN_LENGTH = 50
    MIN_PIDENT = 90.0
    
    print(f"Extracting high-quality coordinates for: {organism}")
    print(f"Filters: length >= {MIN_LENGTH}bp, pident >= {MIN_PIDENT}%")
    print()
    
    # Store coordinates by reference accession
    coords_by_ref = defaultdict(list)
    total_hits = 0
    high_quality = 0
    
    colnames = ['query', 'target', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 
                'taxid', 'org']
    
    with open(input_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 14:
                continue
            
            row_org = fields[13]
            if row_org != organism:
                continue
            
            total_hits += 1
            
            # Parse fields
            target = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            sstart = int(fields[8])
            send = int(fields[9])
            
            # Filter for high quality
            if length >= MIN_LENGTH and pident >= MIN_PIDENT:
                high_quality += 1
                # Normalize coordinates (start < end)
                start, end = (sstart, send) if sstart < send else (send, sstart)
                coords_by_ref[target].append({
                    'start': start,
                    'end': end,
                    'pident': pident,
                    'length': length
                })
    
    print(f"Total hits for {organism}: {total_hits:,}")
    print(f"High-quality hits (≥{MIN_LENGTH}bp, ≥{MIN_PIDENT}%): {high_quality:,} ({100*high_quality/total_hits:.1f}%)")
    print()
    
    if high_quality == 0:
        print("No high-quality hits found.")
        return
    
    # Analyze and output
    print(f"{'Reference Accession':<20} {'# Hits':>8} {'Coordinate Range':>25} {'Span (bp)':>12}")
    print("-" * 80)
    
    output_lines = []
    output_lines.append("accession\tstart\tend\tpident\tlength\n")
    
    for ref_acc in sorted(coords_by_ref.keys()):
        coords = coords_by_ref[ref_acc]
        
        # Get range
        all_starts = [c['start'] for c in coords]
        all_ends = [c['end'] for c in coords]
        min_pos = min(all_starts)
        max_pos = max(all_ends)
        span = max_pos - min_pos
        
        print(f"{ref_acc:<20} {len(coords):>8,} {min_pos:>10,} - {max_pos:>10,} {span:>12,}")
        
        # Sort coordinates by position
        coords_sorted = sorted(coords, key=lambda x: x['start'])
        
        # Write to output
        for c in coords_sorted:
            output_lines.append(f"{ref_acc}\t{c['start']}\t{c['end']}\t{c['pident']:.1f}\t{c['length']}\n")
    
    print()
    
    # Look for clustering
    print("COORDINATE CLUSTERING ANALYSIS:")
    print("-" * 80)
    for ref_acc in sorted(coords_by_ref.keys()):
        coords = sorted(coords_by_ref[ref_acc], key=lambda x: x['start'])
        
        if len(coords) < 2:
            continue
        
        # Calculate gaps between consecutive hits
        gaps = []
        for i in range(len(coords) - 1):
            gap = coords[i+1]['start'] - coords[i]['end']
            gaps.append(gap)
        
        if gaps:
            avg_gap = sum(gaps) / len(gaps)
            max_gap = max(gaps)
            min_gap = min(gaps)
            
            print(f"\n{ref_acc}:")
            print(f"  {len(coords)} hits spanning {coords[0]['start']:,} - {coords[-1]['end']:,}")
            print(f"  Gap between hits: avg={avg_gap:,.0f}bp, min={min_gap:,}bp, max={max_gap:,}bp")
            
            # Check if clustered (small gaps suggest real signal)
            if avg_gap < 1000:
                print(f"  ✓ CLUSTERED: Small gaps suggest contiguous coverage (REAL)")
            elif avg_gap > 100000:
                print(f"  ⚠ SCATTERED: Large gaps suggest random matches (ARTIFACT?)")
    
    # Write output file
    if output_file:
        with open(output_file, 'w') as f:
            f.writelines(output_lines)
        print(f"\nCoordinates written to: {output_file}")
    
    print("\n" + "="*80)
    print("INTERPRETATION:")
    print("="*80)
    print("CLUSTERED hits (small gaps) → likely REAL biological signal")
    print("SCATTERED hits (large gaps) → likely ARTIFACTS (random short matches)")
    print("="*80)

if __name__ == '__main__':
    main()

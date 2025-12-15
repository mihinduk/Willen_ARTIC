#!/usr/bin/env python3
"""
Analyze alignment characteristics for different organisms to detect artifacts.

Real hits should have: high identity + long alignments + good e-values
Artifacts often have: high identity + SHORT alignments (spurious matches)

Usage:
    python3 analyze_organism_hits.py taxonomy_file.txt organism1 organism2 ...
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
# import seaborn as sns

def analyze_organism(df, organism_name):
    """Analyze alignment characteristics for a specific organism."""
    organism_df = df[df['organism'] == organism_name].copy()
    
    if len(organism_df) == 0:
        return None
    
    stats = {
        'organism': organism_name,
        'total_hits': len(organism_df),
        'unique_queries': organism_df['query'].nunique(),
        'mean_pident': organism_df['pident'].mean(),
        'median_pident': organism_df['pident'].median(),
        'mean_length': organism_df['length'].mean(),
        'median_length': organism_df['length'].median(),
        'mean_evalue': organism_df['evalue'].mean(),
        'median_evalue': organism_df['evalue'].median(),
        'mean_bitscore': organism_df['bitscore'].mean(),
        'median_bitscore': organism_df['bitscore'].median(),
        'long_alignments': len(organism_df[organism_df['length'] >= 50]),
        'high_identity': len(organism_df[organism_df['pident'] >= 90]),
        'both_good': len(organism_df[(organism_df['length'] >= 50) & (organism_df['pident'] >= 90)]),
    }
    
    return stats

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 analyze_organism_hits.py taxonomy_file.txt [organism1] [organism2] ...")
        print("\nExample:")
        print("  python3 analyze_organism_hits.py results.txt \"Cyprinus carpio\" \"Severe acute respiratory syndrome coronavirus 2\"")
        sys.exit(1)
    
    input_file = sys.argv[1]
    target_organisms = sys.argv[2:] if len(sys.argv) > 2 else None
    
    print(f"Reading {input_file}...")
    
    # Read the file
    colnames = ['query', 'target', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 
                'taxid', 'organism']
    
    df = pd.read_csv(input_file, sep='\t', names=colnames, dtype={'pident': float, 'length': int, 'evalue': float, 'bitscore': float})
    
    print(f"Loaded {len(df):,} hits")
    print(f"\nTop 10 organisms by hit count:")
    top_organisms = df['organism'].value_counts().head(10)
    for org, count in top_organisms.items():
        pct = 100 * count / len(df)
        print(f"  {org}: {count:,} ({pct:.1f}%)")
    
    # If no specific organisms requested, use top 5
    if not target_organisms:
        target_organisms = top_organisms.head(5).index.tolist()
        print(f"\nAnalyzing top 5 organisms...")
    
    # Analyze each organism
    print(f"\n{'='*80}")
    print("DETAILED ORGANISM ANALYSIS")
    print('='*80)
    
    results = []
    for organism in target_organisms:
        stats = analyze_organism(df, organism)
        if stats:
            results.append(stats)
    
    # Print comparison table
    if results:
        results_df = pd.DataFrame(results)
        
        print("\nAlignment Quality Comparison:")
        print("-" * 80)
        print(f"{'Organism':<50} {'Hits':>10} {'Queries':>10}")
        print("-" * 80)
        for _, row in results_df.iterrows():
            print(f"{row['organism'][:50]:<50} {row['total_hits']:>10,} {row['unique_queries']:>10,}")
        
        print(f"\n{'Organism':<50} {'Med %ID':>8} {'Med Len':>8} {'Med E-val':>12}")
        print("-" * 80)
        for _, row in results_df.iterrows():
            print(f"{row['organism'][:50]:<50} {row['median_pident']:>8.1f} {row['median_length']:>8.0f} {row['median_evalue']:>12.2e}")
        
        print(f"\n{'Organism':<50} {'Long (≥50)':>12} {'High ID (≥90%)':>16} {'Both':>10}")
        print("-" * 80)
        for _, row in results_df.iterrows():
            long_pct = 100 * row['long_alignments'] / row['total_hits']
            high_pct = 100 * row['high_identity'] / row['total_hits']
            both_pct = 100 * row['both_good'] / row['total_hits']
            print(f"{row['organism'][:50]:<50} {row['long_alignments']:>7,} ({long_pct:>4.1f}%) "
                  f"{row['high_identity']:>7,} ({high_pct:>4.1f}%) {row['both_good']:>5,} ({both_pct:>4.1f}%)")
        
        print("\n" + "="*80)
        print("INTERPRETATION:")
        print("="*80)
        print("Real hits should have:")
        print("  - High median alignment length (≥50bp)")
        print("  - High median percent identity (≥90%)")  
        print("  - Good e-values (low)")
        print("  - High % with 'Both' (long + high identity)")
        print("\nArtifacts often have:")
        print("  - SHORT median alignment length (<40bp)")
        print("  - High percent identity but on short sequences")
        print("  - Low % with 'Both' criteria")
        print("="*80)
        
        # Flag suspicious organisms
        print("\n⚠️  POTENTIAL ARTIFACTS:")
        for _, row in results_df.iterrows():
            both_pct = 100 * row['both_good'] / row['total_hits']
            if row['median_length'] < 40 or both_pct < 10:
                print(f"  {row['organism']}: median length={row['median_length']:.0f}bp, "
                      f"only {both_pct:.1f}% have both long+high identity")
        
        print("\n✓  LIKELY REAL:")
        for _, row in results_df.iterrows():
            both_pct = 100 * row['both_good'] / row['total_hits']
            if row['median_length'] >= 50 and both_pct >= 50:
                print(f"  {row['organism']}: median length={row['median_length']:.0f}bp, "
                      f"{both_pct:.1f}% have both long+high identity")

if __name__ == '__main__':
    main()

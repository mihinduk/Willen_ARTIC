#!/usr/bin/env python3
"""
Generate histograms comparing alignment lengths for COVID vs carp (or other organisms)

Usage:
    python3 plot_alignment_length_histograms.py covid_coords.tsv carp_coords.tsv output_prefix
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def load_coords(filename):
    """Load coordinates TSV file."""
    # Columns: accession, start, end, pident, length
    df = pd.read_csv(filename, sep='\t', names=['accession', 'start', 'end', 'pident', 'length'])
    return df

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 plot_alignment_length_histograms.py covid_coords.tsv carp_coords.tsv output_prefix")
        print("\nExample:")
        print("  python3 plot_alignment_length_histograms.py \\")
        print("      Wilen_6_mmseqs_COVID_HQ_coordinates.tsv \\")
        print("      Wilen_6_mmseqs_carp_HQ_coordinates.tsv \\")
        print("      Wilen_6_length_comparison")
        sys.exit(1)
    
    covid_file = sys.argv[1]
    carp_file = sys.argv[2]
    output_prefix = sys.argv[3]
    
    print("="*70)
    print("Alignment Length Comparison")
    print("="*70)
    
    # Load data
    print(f"\nLoading data...")
    print(f"  COVID: {covid_file}")
    covid_df = load_coords(covid_file)
    print(f"    {len(covid_df):,} alignments")
    
    print(f"  Carp: {carp_file}")
    carp_df = load_coords(carp_file)
    print(f"    {len(carp_df):,} alignments")
    
    # Get lengths
    covid_lengths = covid_df['length'].values
    carp_lengths = carp_df['length'].values
    
    # Statistics
    print(f"\n{'='*70}")
    print("STATISTICS")
    print("="*70)
    print(f"\n{'Metric':<30} {'COVID':>15} {'Carp':>15}")
    print("-"*70)
    print(f"{'Number of alignments':<30} {len(covid_lengths):>15,} {len(carp_lengths):>15,}")
    print(f"{'Mean length (bp)':<30} {covid_lengths.mean():>15.1f} {carp_lengths.mean():>15.1f}")
    print(f"{'Median length (bp)':<30} {np.median(covid_lengths):>15.1f} {np.median(carp_lengths):>15.1f}")
    print(f"{'Min length (bp)':<30} {covid_lengths.min():>15.0f} {carp_lengths.min():>15.0f}")
    print(f"{'Max length (bp)':<30} {covid_lengths.max():>15.0f} {carp_lengths.max():>15.0f}")
    print(f"{'Std dev (bp)':<30} {covid_lengths.std():>15.1f} {carp_lengths.std():>15.1f}")
    
    # Percentiles
    print(f"\n{'Percentile':<30} {'COVID':>15} {'Carp':>15}")
    print("-"*70)
    for p in [25, 50, 75, 90, 95, 99]:
        covid_p = np.percentile(covid_lengths, p)
        carp_p = np.percentile(carp_lengths, p)
        print(f"{f'{p}th percentile (bp)':<30} {covid_p:>15.1f} {carp_p:>15.1f}")
    
    # Create figure with two subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Determine shared bin range
    max_len = max(covid_lengths.max(), carp_lengths.max())
    bins = np.arange(50, min(max_len + 10, 300), 5)  # 5bp bins, cap at 300bp for visibility
    
    # Plot 1: COVID histogram
    ax1 = axes[0, 0]
    ax1.hist(covid_lengths, bins=bins, color='#3498db', alpha=0.7, edgecolor='black')
    ax1.axvline(np.median(covid_lengths), color='red', linestyle='--', linewidth=2, 
                label=f'Median: {np.median(covid_lengths):.0f}bp')
    ax1.set_xlabel('Alignment Length (bp)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax1.set_title('COVID-19 Alignment Lengths', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # Plot 2: Carp histogram
    ax2 = axes[0, 1]
    ax2.hist(carp_lengths, bins=bins, color='#e74c3c', alpha=0.7, edgecolor='black')
    ax2.axvline(np.median(carp_lengths), color='blue', linestyle='--', linewidth=2,
                label=f'Median: {np.median(carp_lengths):.0f}bp')
    ax2.set_xlabel('Alignment Length (bp)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax2.set_title('Carp Alignment Lengths', fontsize=12, fontweight='bold')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    
    # Plot 3: Overlaid histograms (normalized)
    ax3 = axes[1, 0]
    ax3.hist(covid_lengths, bins=bins, color='#3498db', alpha=0.5, 
             label=f'COVID (n={len(covid_lengths):,})', density=True, edgecolor='black')
    ax3.hist(carp_lengths, bins=bins, color='#e74c3c', alpha=0.5,
             label=f'Carp (n={len(carp_lengths):,})', density=True, edgecolor='black')
    ax3.set_xlabel('Alignment Length (bp)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Density', fontsize=11, fontweight='bold')
    ax3.set_title('Overlaid Comparison (Normalized)', fontsize=12, fontweight='bold')
    ax3.legend()
    ax3.grid(axis='y', alpha=0.3)
    
    # Plot 4: Cumulative distribution
    ax4 = axes[1, 1]
    covid_sorted = np.sort(covid_lengths)
    carp_sorted = np.sort(carp_lengths)
    covid_cdf = np.arange(1, len(covid_sorted)+1) / len(covid_sorted)
    carp_cdf = np.arange(1, len(carp_sorted)+1) / len(carp_sorted)
    
    ax4.plot(covid_sorted, covid_cdf, color='#3498db', linewidth=2, label='COVID')
    ax4.plot(carp_sorted, carp_cdf, color='#e74c3c', linewidth=2, label='Carp')
    ax4.axhline(0.5, color='gray', linestyle=':', alpha=0.5)
    ax4.axvline(np.median(covid_lengths), color='#3498db', linestyle='--', alpha=0.5)
    ax4.axvline(np.median(carp_lengths), color='#e74c3c', linestyle='--', alpha=0.5)
    ax4.set_xlabel('Alignment Length (bp)', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Cumulative Probability', fontsize=11, fontweight='bold')
    ax4.set_title('Cumulative Distribution Function', fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.grid(alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_png = f"{output_prefix}.png"
    output_pdf = f"{output_prefix}.pdf"
    
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')
    
    print(f"\n{'='*70}")
    print("PLOTS SAVED")
    print("="*70)
    print(f"  {output_png}")
    print(f"  {output_pdf}")
    
    # Interpretation
    print(f"\n{'='*70}")
    print("INTERPRETATION")
    print("="*70)
    
    median_diff = abs(np.median(covid_lengths) - np.median(carp_lengths))
    
    if np.median(covid_lengths) > np.median(carp_lengths) + 10:
        print("\n✓ COVID has significantly LONGER alignments than carp")
        print(f"  → COVID median: {np.median(covid_lengths):.0f}bp")
        print(f"  → Carp median: {np.median(carp_lengths):.0f}bp")
        print(f"  → Difference: {median_diff:.0f}bp")
        print("  → This suggests COVID hits are REAL, carp hits are likely ARTIFACTS")
    elif np.median(carp_lengths) > np.median(covid_lengths) + 10:
        print("\n⚠ Carp has LONGER alignments than COVID (unexpected!)")
        print(f"  → COVID median: {np.median(covid_lengths):.0f}bp")
        print(f"  → Carp median: {np.median(carp_lengths):.0f}bp")
    else:
        print("\n⚠ Similar alignment lengths")
        print(f"  → Both have median ~{np.median(covid_lengths):.0f}bp")
        print("  → Need additional criteria to distinguish real from artifact")
    
    print("\nNote: Even with similar lengths, check the coordinate clustering analysis.")
    print("Real hits should show CLUSTERED coordinates, artifacts are SCATTERED.")
    print("="*70)

if __name__ == '__main__':
    main()

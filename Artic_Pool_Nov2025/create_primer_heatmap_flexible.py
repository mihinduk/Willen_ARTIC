#!/usr/bin/env python3
"""
Create heatmap visualization of missing ARTIC primers across samples
Flexible version that takes filename as argument and handles transposed format

Usage: python3 create_primer_heatmap_flexible.py input_file.tsv [output_prefix]
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def main():
    # Parse command line arguments
    if len(sys.argv) < 2:
        print("Usage: python3 create_primer_heatmap_flexible.py input_file.tsv [output_prefix]")
        print("")
        print("Arguments:")
        print("  input_file.tsv   - Input TSV file with primer data")
        print("  output_prefix    - Optional prefix for output files (default: primer_heatmap)")
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else "primer_heatmap"

    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t', index_col=0)

    # Remove empty columns (trailing tabs)
    df = df.dropna(axis=1, how='all')

    # Remove empty rows
    df = df.dropna(axis=0, how='all')

    print(f"Data shape: {df.shape}")
    print(f"Samples (rows): {list(df.index)}")
    print(f"Number of primers (columns): {len(df.columns)}")

    # Use original orientation: samples on y-axis, primers on x-axis
    df_plot = df

    print(f"\nMissing primers per sample (0 = missing):")
    print((df == 0).sum(axis=1))

    # Determine figure size based on data dimensions
    # samples = rows (height), primers = columns (width)
    n_samples = len(df_plot)
    n_primers = len(df_plot.columns)

    # Calculate height: roughly 0.5 inches per sample, min 6, max 15
    fig_height = max(6, min(15, n_samples * 0.5))
    # Calculate width: roughly 0.05 inches per primer, min 15, max 50
    fig_width = max(15, min(50, n_primers * 0.05))

    print(f"\nFigure size: {fig_width} x {fig_height} inches")

    # Create the heatmap
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Create heatmap with custom colors
    # 0 = missing (orange), 1 = present (blue) - colorblind-friendly palette
    sns.heatmap(df_plot,
                cmap=['#E69F00', '#56B4E9'],  # orange for 0, blue for 1
                cbar_kws={'label': 'Present (1) / Missing (0)', 'ticks': [0.25, 0.75]},
                linewidths=0.5,
                linecolor='lightgray',
                square=False,
                ax=ax)

    # Customize colorbar labels
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticklabels(['Missing', 'Present'])

    # Set labels and title
    ax.set_xlabel('Primer', fontsize=12, fontweight='bold')
    ax.set_ylabel('Sample', fontsize=12, fontweight='bold')
    ax.set_title('Missing ARTIC Primers Across Samples', fontsize=14, fontweight='bold', pad=20)

    # Rotate x-axis labels for better readability (primer names)
    # For many primers, use smaller font and vertical rotation
    if n_primers > 100:
        plt.xticks(rotation=90, ha='center', fontsize=4)
    elif n_primers > 50:
        plt.xticks(rotation=90, ha='center', fontsize=6)
    else:
        plt.xticks(rotation=45, ha='right', fontsize=8)

    # Sample labels on y-axis - easy to read
    plt.yticks(rotation=0, fontsize=10)

    # Tight layout
    plt.tight_layout()

    # Save the figure
    output_png = f"{output_prefix}.png"
    output_pdf = f"{output_prefix}.pdf"

    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')

    print(f"\nHeatmap saved to:")
    print(f"  - {output_png}")
    print(f"  - {output_pdf}")

    # Create summary statistics
    # Remember: 0 = missing, 1 = present
    print("\n" + "="*50)
    print("SUMMARY STATISTICS")
    print("="*50)

    missing_per_sample = (df == 0).sum(axis=1)
    print(f"\nTotal missing primer instances: {(df == 0).sum().sum()}")
    print(f"Average missing primers per sample: {missing_per_sample.mean():.1f}")
    print(f"Sample with most missing primers: {missing_per_sample.idxmax()} ({missing_per_sample.max()} primers)")
    print(f"Sample with fewest missing primers: {missing_per_sample.idxmin()} ({missing_per_sample.min()} primers)")

    print(f"\nPrimers missing in all samples: {(df.sum(axis=0) == 0).sum()}")
    print(f"Primers present in all samples: {(df.sum(axis=0) == len(df)).sum()}")

    # List primers missing in all samples
    all_missing = df.columns[df.sum(axis=0) == 0]
    if len(all_missing) > 0:
        print(f"\nPrimers missing in ALL samples:")
        for primer in all_missing:
            print(f"  - {primer}")

    # List primers present in all samples
    all_present = df.columns[df.sum(axis=0) == len(df)]
    if len(all_present) > 0:
        print(f"\nPrimers present in ALL samples ({len(all_present)} total):")
        if len(all_present) <= 20:
            for primer in all_present:
                print(f"  - {primer}")
        else:
            print(f"  (too many to list - {len(all_present)} primers)")

    # List primers missing most frequently (lowest sums = most 0s)
    missing_counts = (df == 0).sum(axis=0)
    most_missing = missing_counts.sort_values(ascending=False).head(10)
    print(f"\nTop 10 most frequently missing primers:")
    for primer, count in most_missing.items():
        print(f"  - {primer}: {int(count)}/{len(df)} samples")

    print("\nDone!")

if __name__ == "__main__":
    main()

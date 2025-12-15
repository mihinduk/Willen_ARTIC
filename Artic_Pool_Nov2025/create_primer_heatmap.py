#!/usr/bin/env python3
"""
Create heatmap visualization of missing ARTIC primers across samples
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read the TSV file
input_file = "primer_heatmap_table.tsv"
df = pd.read_csv(input_file, sep='\t', index_col=0)

# Remove empty columns (trailing tabs)
df = df.dropna(axis=1, how='all')

# Remove empty rows
df = df.dropna(axis=0, how='all')

print(f"Data shape: {df.shape}")
print(f"Samples: {list(df.columns)}")
print(f"Number of primers: {len(df)}")
print(f"\nMissing primers per sample:")
print(df.sum(axis=0))

# Create the heatmap
fig, ax = plt.subplots(figsize=(10, 14))

# Create heatmap with custom colors
# 0 = present (blue), 1 = missing (orange) - colorblind-friendly palette
sns.heatmap(df,
            cmap=['#56B4E9', '#E69F00'],  # blue for 0, orange for 1
            cbar_kws={'label': 'Missing (1) / Present (0)', 'ticks': [0.25, 0.75]},
            linewidths=0.5,
            linecolor='lightgray',
            square=False,
            ax=ax)

# Customize colorbar labels
colorbar = ax.collections[0].colorbar
colorbar.set_ticklabels(['Present', 'Missing'])

# Set labels and title
ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
ax.set_ylabel('Primer', fontsize=12, fontweight='bold')
ax.set_title('Missing ARTIC Primers Across Samples', fontsize=14, fontweight='bold', pad=20)

# Rotate x-axis labels for better readability
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0, fontsize=8)

# Tight layout
plt.tight_layout()

# Save the figure
output_png = "primer_heatmap.png"
output_pdf = "primer_heatmap.pdf"

plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.savefig(output_pdf, bbox_inches='tight')

print(f"\nHeatmap saved to:")
print(f"  - {output_png}")
print(f"  - {output_pdf}")

# Create summary statistics
print("\n" + "="*50)
print("SUMMARY STATISTICS")
print("="*50)

print(f"\nTotal missing primer instances: {df.sum().sum()}")
print(f"Average missing primers per sample: {df.sum(axis=0).mean():.1f}")
print(f"Samples with most missing primers: {df.sum(axis=0).idxmax()} ({df.sum(axis=0).max()} primers)")
print(f"Samples with fewest missing primers: {df.sum(axis=0).idxmin()} ({df.sum(axis=0).min()} primers)")

print(f"\nPrimers missing in all samples: {(df.sum(axis=1) == len(df.columns)).sum()}")
print(f"Primers missing in no samples: {(df.sum(axis=1) == 0).sum()}")

# List primers missing in all samples
all_missing = df[df.sum(axis=1) == len(df.columns)]
if len(all_missing) > 0:
    print(f"\nPrimers missing in ALL samples:")
    for primer in all_missing.index:
        print(f"  - {primer}")

# List primers missing most frequently
most_missing = df.sum(axis=1).sort_values(ascending=False).head(10)
print(f"\nTop 10 most frequently missing primers:")
for primer, count in most_missing.items():
    print(f"  - {primer}: {int(count)}/{len(df.columns)} samples")

print("\nDone!")

#!/usr/bin/env python3
"""
Generate coverage map from SAM file with gene annotations.

Usage:
    python3 plot_coverage_map.py input.sam

Outputs:
    - input_coverage.png
    - input_coverage.pdf
"""

import sys
import os
import argparse
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import Entrez, SeqIO

# Gene colors for SARS-CoV-2
GENE_COLORS = {
    "orf1ab": "#cd0066",
    "S": "#ff7f50",
    "ORF3a": "#FF8C00",
    "E": "#CD950C",
    "M": "#00FF00",
    "ORF6": "#006400",
    "ORF7a": "#00BFFF",
    "ORF8": "#00008B",
    "N": "#8A2BE2",
    "ORF10": "#000000"
}

# Default email for NCBI Entrez (required)
Entrez.email = "your_email@example.com"


def download_genbank(accession):
    """
    Download GenBank file from NCBI.

    Args:
        accession: GenBank accession (e.g., MN908947.3)

    Returns:
        Path to downloaded GenBank file
    """
    gb_file = f"{accession}.gb"

    if os.path.exists(gb_file):
        print(f"Using existing GenBank file: {gb_file}")
        return gb_file

    print(f"Downloading GenBank file for {accession}...")
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        with open(gb_file, 'w') as f:
            f.write(handle.read())
        handle.close()
        print(f"Downloaded: {gb_file}")
        return gb_file
    except Exception as e:
        print(f"Error downloading GenBank file: {e}")
        return None


def parse_genbank_genes(gb_file):
    """
    Parse gene annotations from GenBank file.

    Args:
        gb_file: Path to GenBank file

    Returns:
        tuple: (genome_length, list of gene dictionaries)
    """
    print(f"Parsing gene annotations from {gb_file}...")

    record = SeqIO.read(gb_file, "genbank")
    genome_length = len(record.seq)

    genes = []
    for feature in record.features:
        if feature.type == "CDS" or feature.type == "mat_peptide":
            # Get gene name
            gene_name = None
            if "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
            elif "product" in feature.qualifiers:
                product = feature.qualifiers["product"][0]
                # Extract gene name from product
                if "ORF" in product.upper():
                    gene_name = product.split()[0]
                elif "spike" in product.lower():
                    gene_name = "S"
                elif "envelope" in product.lower():
                    gene_name = "E"
                elif "membrane" in product.lower():
                    gene_name = "M"
                elif "nucleocapsid" in product.lower():
                    gene_name = "N"

            if gene_name:
                # Normalize gene name
                if gene_name.lower() == "orf1a" or gene_name.lower() == "orf1b":
                    gene_name = "orf1ab"

                # Get start and end positions
                start = int(feature.location.start) + 1  # Convert to 1-based
                end = int(feature.location.end)

                genes.append({
                    'name': gene_name,
                    'start': start,
                    'end': end,
                    'product': feature.qualifiers.get("product", [""])[0]
                })

    # Merge overlapping orf1ab entries
    orf1ab_genes = [g for g in genes if g['name'] == 'orf1ab']
    if len(orf1ab_genes) > 1:
        merged_start = min(g['start'] for g in orf1ab_genes)
        merged_end = max(g['end'] for g in orf1ab_genes)
        genes = [g for g in genes if g['name'] != 'orf1ab']
        genes.append({
            'name': 'orf1ab',
            'start': merged_start,
            'end': merged_end,
            'product': 'ORF1ab polyprotein'
        })

    # Sort by start position
    genes.sort(key=lambda x: x['start'])

    print(f"Found {len(genes)} genes:")
    for gene in genes:
        print(f"  {gene['name']}: {gene['start']}-{gene['end']} ({gene['product']})")

    return genome_length, genes


def parse_cigar(cigar_string):
    """
    Parse CIGAR string to determine reference positions covered.

    Returns list of (operation, length) tuples.
    Operations that consume reference: M, D, N, =, X
    Operations that don't: I, S, H, P
    """
    import re
    return [(op, int(length)) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)]


def parse_sam_coverage(sam_file, min_mapq=30):
    """
    Parse full SAM file and calculate coverage depth at each position.
    Uses CIGAR string to determine actual reference positions covered.

    Expected format (tab-separated, standard SAM):
    QNAME FLAG RNAME POS MAPQ CIGAR SEQ QUAL [optional fields...]

    Args:
        sam_file: Path to SAM file
        min_mapq: Minimum mapping quality (default: 30)

    Returns:
        dict: Position -> depth
    """
    print(f"Parsing SAM file: {sam_file}...")
    print(f"Minimum MAPQ filter: {min_mapq}")
    print(f"Using CIGAR strings to determine actual coverage")

    coverage = defaultdict(int)
    total_reads = 0
    passed_reads = 0
    failed_mapq = 0
    failed_unmapped = 0
    total_bases_covered = 0

    with open(sam_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # Skip header lines
            if line.startswith('@'):
                continue

            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            if len(fields) < 11:
                print(f"Warning: Line {line_num} has fewer than 11 fields (standard SAM), skipping")
                continue

            total_reads += 1

            # Parse standard SAM fields
            try:
                qname = fields[0]
                flag = int(fields[1])
                rname = fields[2]
                pos = int(fields[3])  # 1-based leftmost position
                mapq = int(fields[4])
                cigar = fields[5]
            except (ValueError, IndexError) as e:
                print(f"Warning: Error parsing line {line_num}: {e}")
                continue

            # Skip unmapped reads (flag & 0x4)
            if flag & 0x4:
                failed_unmapped += 1
                continue

            # Filter by mapping quality
            if mapq < min_mapq:
                failed_mapq += 1
                continue

            passed_reads += 1

            # Parse CIGAR to determine covered positions
            cigar_ops = parse_cigar(cigar)
            current_pos = pos

            for op, length in cigar_ops:
                if op in ['M', '=', 'X']:  # Match/mismatch - covers reference
                    for i in range(length):
                        coverage[current_pos + i] += 1
                    current_pos += length
                    total_bases_covered += length
                elif op in ['D', 'N']:  # Deletion/skip - advances reference but no coverage
                    current_pos += length
                elif op in ['I', 'S', 'H', 'P']:  # Insertion/soft-clip/hard-clip/padding - no reference advance
                    pass  # Don't advance reference position

    print(f"\nRead filtering statistics:")
    print(f"  Total reads: {total_reads:,}")
    print(f"  Unmapped reads: {failed_unmapped:,} ({100*failed_unmapped/total_reads:.1f}%)" if total_reads > 0 else "  Unmapped reads: 0")
    print(f"  Failed MAPQ<{min_mapq}: {failed_mapq:,} ({100*failed_mapq/total_reads:.1f}%)" if total_reads > 0 else f"  Failed MAPQ<{min_mapq}: 0")
    print(f"  Passed filters: {passed_reads:,} ({100*passed_reads/total_reads:.1f}%)" if total_reads > 0 else "  Passed filters: 0")
    if passed_reads > 0:
        print(f"  Average mapped bases per read: {total_bases_covered/passed_reads:.1f} bp")

    if coverage:
        covered_positions = len(coverage)
        max_depth = max(coverage.values())
        avg_depth = sum(coverage.values()) / len(coverage)
        min_pos = min(coverage.keys())
        max_pos = max(coverage.keys())
        print(f"\nCoverage statistics:")
        print(f"  Coverage range: {min_pos:,} - {max_pos:,} bp")
        print(f"  Covered positions: {covered_positions:,}")
        print(f"  Max depth: {max_depth:,}x")
        print(f"  Average depth: {avg_depth:.1f}x")
    else:
        print("\nWarning: No coverage data after filtering")

    return coverage


def plot_coverage_map(coverage, genes, genome_length, output_prefix, log_scale=True):
    """
    Create coverage map visualization with gene annotations.
    Shows entire genome including 5' and 3' UTR regions.

    Args:
        coverage: Dict of position -> depth
        genes: List of gene dictionaries
        genome_length: Total genome length
        output_prefix: Output file prefix
        log_scale: If True, plot log10(depth+1), otherwise linear (default: True)
    """
    print("\nCreating coverage visualization...")
    print(f"  Genome length: {genome_length:,} bp")
    print(f"  Scale: {'log10(depth+1)' if log_scale else 'linear'}")

    # Calculate 5' and 3' UTR regions
    first_gene_start = min(g['start'] for g in genes) if genes else 1
    last_gene_end = max(g['end'] for g in genes) if genes else genome_length

    utr5_length = first_gene_start - 1
    utr3_length = genome_length - last_gene_end

    print(f"  5' UTR: 1-{first_gene_start-1} ({utr5_length} bp)")
    print(f"  Coding: {first_gene_start}-{last_gene_end} ({last_gene_end-first_gene_start+1} bp)")
    print(f"  3' UTR: {last_gene_end+1}-{genome_length} ({utr3_length} bp)")

    # Prepare data for plotting - ensure we cover the ENTIRE genome
    positions = list(range(1, genome_length + 1))
    depths = [coverage.get(pos, 0) for pos in positions]

    # Apply log transformation if requested
    if log_scale:
        depths = np.log10(np.array(depths) + 1)

    # Create figure with two subplots
    fig, (ax_cov, ax_genes) = plt.subplots(
        2, 1,
        figsize=(16, 8),
        gridspec_kw={'height_ratios': [3, 1]},
        sharex=True
    )

    # Plot coverage
    ax_cov.fill_between(positions, depths, color='steelblue', alpha=0.7)
    ax_cov.plot(positions, depths, color='navy', linewidth=0.5, alpha=0.5)

    # Add 50X dashed line
    if log_scale:
        depth_50x = np.log10(50 + 1)
        ax_cov.axhline(y=depth_50x, color='red', linestyle='--', linewidth=2, label='50X depth')
    else:
        ax_cov.axhline(y=50, color='red', linestyle='--', linewidth=2, label='50X depth')

    # Format coverage plot
    if log_scale:
        ax_cov.set_ylabel('Read Depth [log10(depth+1)]', fontsize=12, fontweight='bold')
    else:
        ax_cov.set_ylabel('Read Depth', fontsize=12, fontweight='bold')
    ax_cov.set_ylim(bottom=0)
    ax_cov.grid(True, alpha=0.3)
    ax_cov.legend(loc='upper right', fontsize=10)
    ax_cov.set_title('Coverage Map with Gene Annotations', fontsize=14, fontweight='bold', pad=20)

    # Plot genes
    ax_genes.set_ylim(0, 1)
    ax_genes.set_xlim(0, genome_length)

    # Highlight 5' and 3' UTR regions with light gray shading
    if utr5_length > 0:
        utr5_rect = patches.Rectangle(
            (0, 0), first_gene_start - 1, 1,
            linewidth=0, facecolor='lightgray', alpha=0.3
        )
        ax_genes.add_patch(utr5_rect)
        # Label 5' UTR - staggered at bottom
        if utr5_length > 100:  # Only label if visible
            ax_genes.text(
                (first_gene_start - 1) / 2, 0.15, "5' UTR",
                ha='center', va='center',
                fontsize=8, style='italic', color='gray'
            )

    if utr3_length > 0:
        utr3_rect = patches.Rectangle(
            (last_gene_end, 0), genome_length - last_gene_end, 1,
            linewidth=0, facecolor='lightgray', alpha=0.3
        )
        ax_genes.add_patch(utr3_rect)
        # Label 3' UTR - staggered at top
        if utr3_length > 100:  # Only label if visible
            ax_genes.text(
                last_gene_end + utr3_length / 2, 0.85, "3' UTR",
                ha='center', va='center',
                fontsize=8, style='italic', color='gray'
            )

    for i, gene in enumerate(genes):
        gene_name = gene['name']
        start = gene['start']
        end = gene['end']

        # Get color (case-insensitive lookup)
        color = GENE_COLORS.get(gene_name, '#808080')
        if color == '#808080':  # Try uppercase
            color = GENE_COLORS.get(gene_name.upper(), '#808080')

        # Draw gene box
        width = end - start
        rect = patches.Rectangle(
            (start, 0.2), width, 0.6,
            linewidth=1, edgecolor='black', facecolor=color, alpha=0.8
        )
        ax_genes.add_patch(rect)

        # Stagger gene labels to avoid overlap
        # Alternate between three vertical positions
        label_positions = [0.5, 0.3, 0.7]  # center, lower, upper
        label_y = label_positions[i % 3]

        label_pos = (start + end) / 2
        ax_genes.text(
            label_pos, label_y, gene_name,
            ha='center', va='center',
            fontsize=9, fontweight='bold',
            color='white' if gene_name in ['ORF6', 'ORF8', 'N', 'ORF10'] else 'black'
        )

    # Format gene plot
    ax_genes.set_xlabel('Genome Position (bp)', fontsize=12, fontweight='bold')
    ax_genes.set_ylabel('Genes', fontsize=12, fontweight='bold')
    ax_genes.set_yticks([])
    ax_genes.spines['left'].set_visible(False)
    ax_genes.spines['right'].set_visible(False)
    ax_genes.spines['top'].set_visible(False)

    # Format x-axis
    ax_genes.ticklabel_format(style='plain', axis='x')

    # Adjust layout
    plt.tight_layout()

    # Save outputs
    png_file = f"{output_prefix}_coverage.png"
    pdf_file = f"{output_prefix}_coverage.pdf"

    plt.savefig(png_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {png_file}")

    plt.savefig(pdf_file, bbox_inches='tight')
    print(f"Saved: {pdf_file}")

    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Generate coverage map from SAM file with gene annotations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use with SARS-CoV-2 mappings (auto-download annotation, log scale by default)
  python3 plot_coverage_map.py Handley_A2690_54317_Wilen_9_mappings.nonsupplemental.sam --accession MN908947.3

  # Adjust quality filter (default MAPQ=30)
  python3 plot_coverage_map.py mappings.sam --accession MN908947.3 --min-mapq 60

  # Use linear scale instead of log scale
  python3 plot_coverage_map.py mappings.sam --accession MN908947.3 --linear

  # Use existing GenBank file
  python3 plot_coverage_map.py mappings.sam --genbank MN908947.3.gb

Input Format:
  Expects standard SAM format (tab-separated, 11+ columns):
    QNAME FLAG RNAME POS MAPQ CIGAR SEQ QUAL [optional fields...]

  Uses CIGAR string to determine actual reference positions covered.
  Properly handles soft/hard clipping, insertions, deletions, etc.

  Generated from BWA-MEM alignments using:
    bwa mem -t 12 reference.fasta R1.fastq.gz R2.fastq.gz > mappings.sam
    samtools view -h -F 2048 mappings.sorted.bam > mappings.nonsupplemental.sam

Output:
  Creates two files showing coverage across entire genome including 5' and 3' UTRs:
    - {sample}_coverage.png (high resolution, 300 DPI)
    - {sample}_coverage.pdf (vector graphics)

  By default, plots log10(depth+1) to better visualize coverage dynamics.
  Use --linear flag for linear scale.
        """
    )

    parser.add_argument('sam_file', help='Input SAM file (standard format with CIGAR strings)')
    parser.add_argument('--accession',
                       help='GenBank accession for gene annotations (e.g., MN908947.3)')
    parser.add_argument('--genbank', help='Use existing GenBank file instead of downloading')
    parser.add_argument('--email', help='Email for NCBI Entrez (required for downloads)')
    parser.add_argument('--min-mapq', type=int, default=30,
                       help='Minimum mapping quality (default: 30, recommended: 20-60)')
    parser.add_argument('--linear', action='store_true',
                       help='Use linear scale instead of log10(depth+1) (default: log scale)')

    args = parser.parse_args()

    # Check SAM file exists
    if not os.path.exists(args.sam_file):
        print(f"Error: SAM file not found: {args.sam_file}")
        sys.exit(1)

    # Require either accession or genbank
    if not args.accession and not args.genbank:
        print("Error: Must provide either --accession or --genbank")
        print("\nExamples:")
        print("  python3 plot_coverage_map.py sample.sam --accession MN908947.3")
        print("  python3 plot_coverage_map.py sample.sam --genbank MN908947.3.gb --min-mapq 30")
        sys.exit(1)

    # Set email if provided
    if args.email:
        Entrez.email = args.email

    # Get or download GenBank file
    if args.genbank:
        gb_file = args.genbank
        if not os.path.exists(gb_file):
            print(f"Error: GenBank file not found: {gb_file}")
            sys.exit(1)
    else:
        gb_file = download_genbank(args.accession)
        if not gb_file:
            print("Failed to download GenBank file")
            sys.exit(1)

    # Parse gene annotations
    genome_length, genes = parse_genbank_genes(gb_file)

    # Parse SAM file for coverage with quality filtering
    coverage = parse_sam_coverage(args.sam_file, min_mapq=args.min_mapq)

    # Generate output prefix from SAM filename
    output_prefix = os.path.splitext(args.sam_file)[0]

    # Create visualization
    log_scale = not args.linear  # Default is log scale unless --linear is specified
    plot_coverage_map(coverage, genes, genome_length, output_prefix, log_scale=log_scale)

    print("\n" + "="*60)
    print("SUCCESS!")
    print("="*60)
    print(f"Generated coverage maps:")
    print(f"  {output_prefix}_coverage.png")
    print(f"  {output_prefix}_coverage.pdf")
    print("\n" + "="*60)
    print("Veritas numquam perit")
    print("Truth never perishes")
    print("="*60)


if __name__ == '__main__':
    main()

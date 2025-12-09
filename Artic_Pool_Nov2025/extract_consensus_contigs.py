#!/usr/bin/env python3
"""
Extract consensus contigs from SAM file for BLAST validation.
Generates consensus sequence from mapped reads and splits into contigs at coverage gaps.

Usage:
    python3 extract_consensus_contigs.py input.sam reference.fasta output.fasta [--min-depth 10] [--min-contig-length 100]
"""

import sys
import re
import argparse
from collections import defaultdict, Counter

def parse_cigar(cigar_string, start_pos, sequence):
    """Parse CIGAR and return positions with their bases."""
    positions = {}
    current_pos = start_pos
    seq_pos = 0

    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    for length, op in ops:
        length = int(length)
        if op in ['M', '=', 'X']:  # Match/mismatch
            for i in range(length):
                positions[current_pos + i] = sequence[seq_pos + i]
            seq_pos += length
            current_pos += length
        elif op in ['D', 'N']:  # Deletion - advance reference
            current_pos += length
        elif op in ['I', 'S']:  # Insertion/soft-clip - advance sequence
            seq_pos += length

    return positions

def parse_sam_pileup(sam_file, min_mapq=30):
    """Parse SAM and create base pileup at each position."""
    pileup = defaultdict(list)

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            flag = int(fields[1])
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]

            # Skip unmapped or low quality
            if flag & 0x4 or mapq < min_mapq:
                continue

            # Get bases at each position
            positions = parse_cigar(cigar, pos, seq)
            for p, base in positions.items():
                pileup[p].append(base.upper())

    return pileup

def call_consensus(pileup, min_depth=10):
    """Call consensus base at each position with sufficient depth."""
    consensus = {}

    for pos in sorted(pileup.keys()):
        bases = pileup[pos]
        if len(bases) >= min_depth:
            # Call most common base
            base_counts = Counter(bases)
            consensus_base = base_counts.most_common(1)[0][0]
            consensus[pos] = consensus_base

    return consensus

def extract_contigs(consensus, min_contig_length=100):
    """Extract contiguous regions from consensus as separate contigs."""
    if not consensus:
        return []

    contigs = []
    positions = sorted(consensus.keys())

    current_contig_start = positions[0]
    current_contig_seq = [consensus[positions[0]]]

    for i in range(1, len(positions)):
        pos = positions[i]
        prev_pos = positions[i-1]

        if pos == prev_pos + 1:  # Contiguous
            current_contig_seq.append(consensus[pos])
        else:  # Gap detected
            # Save current contig if long enough
            if len(current_contig_seq) >= min_contig_length:
                contigs.append((
                    current_contig_start,
                    prev_pos,
                    ''.join(current_contig_seq)
                ))
            # Start new contig
            current_contig_start = pos
            current_contig_seq = [consensus[pos]]

    # Don't forget the last contig
    if len(current_contig_seq) >= min_contig_length:
        contigs.append((
            current_contig_start,
            positions[-1],
            ''.join(current_contig_seq)
        ))

    return contigs

def main():
    parser = argparse.ArgumentParser(
        description='Extract consensus contigs from SAM file for BLAST validation'
    )
    parser.add_argument('sam_file', help='Input SAM file')
    parser.add_argument('output_fasta', help='Output FASTA file with contigs')
    parser.add_argument('--min-depth', type=int, default=10,
                        help='Minimum depth for consensus calling (default: 10)')
    parser.add_argument('--min-contig-length', type=int, default=100,
                        help='Minimum contig length to output (default: 100bp)')
    parser.add_argument('--min-mapq', type=int, default=30,
                        help='Minimum MAPQ for filtering reads (default: 30)')

    args = parser.parse_args()

    print(f"Processing: {args.sam_file}")
    print(f"Min depth for consensus: {args.min_depth}X")
    print(f"Min contig length: {args.min_contig_length}bp")
    print(f"Min MAPQ: {args.min_mapq}\n")

    # Parse SAM and build pileup
    print("Building base pileup from SAM file...")
    pileup = parse_sam_pileup(args.sam_file, args.min_mapq)
    print(f"Positions with coverage: {len(pileup):,}\n")

    # Call consensus
    print("Calling consensus bases...")
    consensus = call_consensus(pileup, args.min_depth)
    print(f"Consensus positions (≥{args.min_depth}X): {len(consensus):,}\n")

    # Extract contigs
    print("Extracting contigs...")
    contigs = extract_contigs(consensus, args.min_contig_length)
    print(f"Contigs found (≥{args.min_contig_length}bp): {len(contigs)}\n")

    # Write output
    with open(args.output_fasta, 'w') as f:
        for i, (start, end, seq) in enumerate(contigs, 1):
            contig_name = f"contig_{i}_pos{start}-{end}_len{len(seq)}"
            f.write(f">{contig_name}\n")
            # Write sequence in 80-character lines
            for j in range(0, len(seq), 80):
                f.write(f"{seq[j:j+80]}\n")

    print("=" * 70)
    print("CONTIGS EXTRACTED")
    print("=" * 70)
    for i, (start, end, seq) in enumerate(contigs, 1):
        print(f"Contig {i}: positions {start:,}-{end:,} ({len(seq):,} bp)")

    print(f"\nTotal bases in contigs: {sum(len(seq) for _, _, seq in contigs):,}")
    print(f"\nWrote contigs to: {args.output_fasta}")
    print("\nReady for BLAST against nt database!")

if __name__ == '__main__':
    main()

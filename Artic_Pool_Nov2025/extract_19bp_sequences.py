#!/usr/bin/env python3
"""
Extract sequences from SAM file with exactly 19bp mapped length.
Creates two FASTA files:
  1. Full sequences (including soft-clipped bases)
  2. Aligned-only sequences (only the 19bp that mapped)

Usage:
    python3 extract_19bp_sequences.py input.sam output_prefix
"""

import sys
import re

def parse_cigar_length(cigar):
    """Calculate mapped length from CIGAR."""
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    mapped_length = sum(int(num) for num, op in ops if op in ['M', '=', 'X'])
    return mapped_length

def extract_aligned_sequence(seq, cigar):
    """
    Extract only the aligned portion of the sequence based on CIGAR.
    Skips soft-clipped (S), hard-clipped (H), and inserted (I) bases.
    """
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    aligned_seq = []
    seq_pos = 0  # Position in the read sequence

    for length, op in ops:
        length = int(length)

        if op in ['M', '=', 'X']:  # Match/mismatch - part of alignment
            aligned_seq.append(seq[seq_pos:seq_pos + length])
            seq_pos += length
        elif op in ['I', 'S']:  # Insertion or soft-clip - skip these bases
            seq_pos += length
        elif op in ['D', 'N']:  # Deletion/skip - doesn't consume query
            pass
        elif op == 'H':  # Hard clip - already removed from SEQ
            pass

    return ''.join(aligned_seq)

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 extract_19bp_sequences.py input.sam output_prefix")
        print("\nExample:")
        print("  python3 extract_19bp_sequences.py sample.sam sample_19bp")
        print("\nOutputs:")
        print("  sample_19bp_full.fasta       - Full sequences (with soft-clips)")
        print("  sample_19bp_aligned.fasta    - Only aligned 19bp")
        sys.exit(1)

    sam_file = sys.argv[1]
    output_prefix = sys.argv[2]

    output_full = f"{output_prefix}_full.fasta"
    output_aligned = f"{output_prefix}_aligned.fasta"

    print(f"Processing: {sam_file}")
    print(f"Extracting reads with exactly 19bp mapped length\n")

    total_reads = 0
    target_reads = 0

    with open(sam_file, 'r') as infile, \
         open(output_full, 'w') as full_out, \
         open(output_aligned, 'w') as aligned_out:

        for line in infile:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            total_reads += 1

            read_name = fields[0]
            flag = int(fields[1])
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]

            # Skip unmapped
            if flag & 0x4:
                continue

            # Check if mapped length is exactly 19bp
            mapped_length = parse_cigar_length(cigar)

            if mapped_length == 19:
                target_reads += 1

                # Extract aligned portion
                aligned_seq = extract_aligned_sequence(seq, cigar)

                # Write full sequence (including soft-clips)
                full_out.write(f">{read_name}_pos{pos}_mapq{mapq}_full\n")
                full_out.write(f"{seq}\n")

                # Write aligned-only sequence
                aligned_out.write(f">{read_name}_pos{pos}_mapq{mapq}_aligned\n")
                aligned_out.write(f"{aligned_seq}\n")

                # Show first 5 examples
                if target_reads <= 5:
                    print(f"Example {target_reads}:")
                    print(f"  Read: {read_name}")
                    print(f"  Position: {pos}, MAPQ: {mapq}")
                    print(f"  CIGAR: {cigar}")
                    print(f"  Full seq ({len(seq)}bp): {seq[:60]}...")
                    print(f"  Aligned ({len(aligned_seq)}bp): {aligned_seq}")
                    print()

            if total_reads % 100000 == 0:
                print(f"  Processed {total_reads:,} reads, found {target_reads:,} with 19bp mapped...")

    print(f"\n{'='*60}")
    print(f"EXTRACTION COMPLETE")
    print(f"{'='*60}")
    print(f"Total reads processed: {total_reads:,}")
    print(f"Reads with exactly 19bp mapped: {target_reads:,} ({100*target_reads/total_reads:.2f}%)")
    print(f"\nOutput files:")
    print(f"  {output_full} - Full sequences ({target_reads:,} reads)")
    print(f"  {output_aligned} - Aligned-only sequences ({target_reads:,} reads)")

if __name__ == '__main__':
    main()

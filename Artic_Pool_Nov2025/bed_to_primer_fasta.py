#!/usr/bin/env python3
"""
Extract primer sequences from genome FASTA using BED file coordinates.

Usage:
    python3 bed_to_primer_fasta.py primers.bed genome.fasta > primers.fasta
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq

def parse_bed(bed_file):
    """Parse BED file and return list of primer records."""
    primers = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 6:
                # Basic BED format: chr, start, end, name
                if len(fields) >= 4:
                    chrom = fields[0]
                    start = int(fields[1])  # 0-based
                    end = int(fields[2])
                    name = fields[3]
                    strand = '+'  # Default to plus if not specified
                else:
                    continue
            else:
                # Extended BED: chr, start, end, name, score, strand
                chrom = fields[0]
                start = int(fields[1])  # 0-based
                end = int(fields[2])
                name = fields[3]
                strand = fields[5]

            primers.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': name,
                'strand': strand
            })

    return primers

def extract_sequences(primers, genome_file):
    """Extract primer sequences from genome."""
    # Load genome
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta'))

    primer_seqs = []

    for primer in primers:
        chrom = primer['chrom']

        # Try to find matching chromosome
        if chrom in genome:
            seq_record = genome[chrom]
        else:
            # Try without version number (e.g., MN908947 instead of MN908947.3)
            chrom_base = chrom.split('.')[0]
            matching = [k for k in genome.keys() if k.startswith(chrom_base)]
            if matching:
                seq_record = genome[matching[0]]
            else:
                print(f"Warning: Chromosome {chrom} not found in genome", file=sys.stderr)
                continue

        # Extract sequence (BED is 0-based, half-open)
        start = primer['start']
        end = primer['end']
        seq = seq_record.seq[start:end]

        # Reverse complement if on minus strand
        if primer['strand'] == '-':
            seq = seq.reverse_complement()

        primer_seqs.append({
            'name': primer['name'],
            'seq': str(seq),
            'length': len(seq),
            'strand': primer['strand'],
            'coords': f"{chrom}:{start}-{end}"
        })

    return primer_seqs

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 bed_to_primer_fasta.py primers.bed genome.fasta > primers.fasta", file=sys.stderr)
        sys.exit(1)

    bed_file = sys.argv[1]
    genome_file = sys.argv[2]

    print(f"Reading BED file: {bed_file}", file=sys.stderr)
    primers = parse_bed(bed_file)
    print(f"Found {len(primers)} primers", file=sys.stderr)

    print(f"Extracting sequences from: {genome_file}", file=sys.stderr)
    primer_seqs = extract_sequences(primers, genome_file)
    print(f"Extracted {len(primer_seqs)} primer sequences", file=sys.stderr)

    # Output FASTA
    for primer in primer_seqs:
        print(f">{primer['name']} {primer['coords']} strand={primer['strand']} length={primer['length']}")
        print(primer['seq'])

    print(f"\nWrote {len(primer_seqs)} primers to stdout", file=sys.stderr)

if __name__ == '__main__':
    main()

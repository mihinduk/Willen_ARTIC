#!/usr/bin/env python3
"""
Parse BLAST results and classify hits.
- Non-COVID: classify by taxonomy (human/bacterial/other), create FASTAs
- COVID hits: map to genes and create TSV

Usage:
    python3 parse_blast_results.py blast_results.txt reads.fasta sam_metadata.tsv output_prefix
"""

import sys
import re
from collections import defaultdict, Counter

# SARS-CoV-2 MN908947.3 gene annotations (1-based coordinates)
SARSCOV2_GENES = [
    {'name': 'orf1ab', 'start': 266, 'end': 21555, 'product': 'ORF1ab polyprotein'},
    {'name': 'S', 'start': 21563, 'end': 25384, 'product': 'surface glycoprotein'},
    {'name': 'ORF3a', 'start': 25393, 'end': 26220, 'product': 'ORF3a protein'},
    {'name': 'E', 'start': 26245, 'end': 26472, 'product': 'envelope protein'},
    {'name': 'M', 'start': 26523, 'end': 27191, 'product': 'membrane glycoprotein'},
    {'name': 'ORF6', 'start': 27202, 'end': 27387, 'product': 'ORF6 protein'},
    {'name': 'ORF7a', 'start': 27394, 'end': 27759, 'product': 'ORF7a protein'},
    {'name': 'ORF7b', 'start': 27756, 'end': 27887, 'product': 'ORF7b protein'},
    {'name': 'ORF8', 'start': 27894, 'end': 28259, 'product': 'ORF8 protein'},
    {'name': 'N', 'start': 28274, 'end': 29533, 'product': 'nucleocapsid phosphoprotein'},
    {'name': 'ORF10', 'start': 29558, 'end': 29674, 'product': 'ORF10 protein'},
]

def map_position_to_genes(start, end):
    """Map read positions to overlapping genes."""
    genes_hit = []
    for gene in SARSCOV2_GENES:
        # Check if read overlaps with gene
        if not (end < gene['start'] or start > gene['end']):
            genes_hit.append(gene['name'])

    if not genes_hit:
        # Check if in UTRs
        if end < SARSCOV2_GENES[0]['start']:
            genes_hit.append("5'_UTR")
        elif start > SARSCOV2_GENES[-1]['end']:
            genes_hit.append("3'_UTR")
        else:
            genes_hit.append("intergenic")

    return ",".join(genes_hit)

def parse_read_id(read_id):
    """Extract mapping coordinates from read ID."""
    # Format: readname_pos123-456_len78_mapq30
    match = re.search(r'pos(\d+)-(\d+)_len(\d+)', read_id)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        length = int(match.group(3))
        return start, end, length
    return None, None, None

def load_fasta_sequences(fasta_file):
    """Load FASTA sequences into dictionary."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences

def classify_taxonomy(sciname):
    """Classify organism by scientific name."""
    sciname_lower = sciname.lower()

    # COVID variants
    if 'sars' in sciname_lower and 'cov' in sciname_lower:
        return 'covid'
    if 'severe acute respiratory syndrome' in sciname_lower:
        return 'covid'

    # Human
    if 'homo sapiens' in sciname_lower or 'human' in sciname_lower:
        return 'human'

    # Bacterial
    bacterial_terms = ['bacteria', 'bacillus', 'escherichia', 'staphylococcus',
                      'streptococcus', 'pseudomonas', 'mycobacterium', 'clostridium',
                      'salmonella', 'shigella', 'vibrio', 'listeria']
    if any(term in sciname_lower for term in bacterial_terms):
        return 'bacterial'

    # Viral (non-COVID)
    if 'virus' in sciname_lower or 'phage' in sciname_lower:
        return 'viral'

    # Fungal
    if 'fungus' in sciname_lower or 'fungi' in sciname_lower or 'candida' in sciname_lower:
        return 'fungal'

    return 'other'

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 parse_blast_results.py blast_results.txt reads.fasta sam_metadata.tsv output_prefix")
        print("\nBLAST results should be in format 6 with columns:")
        print("  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames")
        print("\nExample BLAST command:")
        print("  blastn -query reads.fasta -db nt -num_threads 16 \\")
        print("    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \\")
        print("    -max_target_seqs 5 -out blast_results.txt")
        sys.exit(1)

    blast_file = sys.argv[1]
    fasta_file = sys.argv[2]
    metadata_file = sys.argv[3]
    output_prefix = sys.argv[4]

    print(f"Processing BLAST results: {blast_file}")
    print(f"Loading sequences from: {fasta_file}")
    print(f"Loading metadata from: {metadata_file}\n")

    # Load sequences
    sequences = load_fasta_sequences(fasta_file)
    print(f"Loaded {len(sequences):,} sequences\n")

    # Parse BLAST results - collect all hits per read
    read_hits = defaultdict(list)  # qseqid -> list of hits

    with open(blast_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 14:
                continue

            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            evalue = float(fields[10])
            bitscore = float(fields[11])
            staxids = fields[12] if len(fields) > 12 else "N/A"
            sscinames = fields[13] if len(fields) > 13 else "Unknown"

            read_hits[qseqid].append({
                'sseqid': sseqid,
                'pident': pident,
                'evalue': evalue,
                'bitscore': bitscore,
                'sscinames': sscinames,
                'taxonomy': classify_taxonomy(sscinames)
            })

    print(f"Loaded BLAST results for {len(read_hits):,} reads\n")

    # Analyze hits and classify reads
    covid_hits = []
    non_covid_hits = defaultdict(list)
    taxonomy_counts = Counter()
    ambiguous_reads = []
    top5_agreement = Counter()  # Track how often top 5 hits agree

    for qseqid, hits in read_hits.items():
        if not hits:
            continue

        # Sort by bitscore (best first) - should already be sorted but ensure
        hits = sorted(hits, key=lambda x: x['bitscore'], reverse=True)

        # Use top hit for classification
        top_hit = hits[0]
        taxonomy = top_hit['taxonomy']
        taxonomy_counts[taxonomy] += 1

        # Check agreement among top hits
        top_taxonomies = [h['taxonomy'] for h in hits[:5]]
        unique_taxonomies = set(top_taxonomies)

        if len(unique_taxonomies) == 1:
            top5_agreement['unanimous'] += 1
        elif len(unique_taxonomies) == 2:
            top5_agreement['mostly_agree'] += 1
        else:
            top5_agreement['ambiguous'] += 1
            ambiguous_reads.append((qseqid, top_taxonomies, [h['sscinames'] for h in hits[:5]]))

        # Process based on top hit taxonomy
        if taxonomy == 'covid':
            # Extract coordinates from read ID
            start, end, length = parse_read_id(qseqid)
            if start and end:
                genes = map_position_to_genes(start, end)
                # Extract original read name (before _pos...)
                read_name = re.sub(r'_pos\d+-\d+_len\d+_mapq\d+', '', qseqid)

                # Check if all top 5 are COVID
                all_covid = all(h['taxonomy'] == 'covid' for h in hits[:min(5, len(hits))])

                covid_hits.append({
                    'read_name': read_name,
                    'mapped_length': length,
                    'start': start,
                    'end': end,
                    'genes': genes,
                    'pident': top_hit['pident'],
                    'evalue': top_hit['evalue'],
                    'top5_all_covid': all_covid,
                    'num_hits': len(hits)
                })
        else:
            non_covid_hits[taxonomy].append((qseqid, top_hit['sscinames'],
                                             top_hit['pident'], top_hit['evalue']))

    print(f"Processed {len(read_hits):,} reads with BLAST hits\n")

    # Output summary
    print("="*70)
    print("TAXONOMY SUMMARY (based on top hit)")
    print("="*70)
    total_reads = len(read_hits)
    for taxonomy in sorted(taxonomy_counts.keys()):
        count = taxonomy_counts[taxonomy]
        pct = 100 * count / total_reads if total_reads else 0
        print(f"{taxonomy.upper()}: {count:,} reads ({pct:.2f}%)")
    print()

    # Top-5 agreement summary
    print("="*70)
    print("TOP-5 HIT AGREEMENT")
    print("="*70)
    for agreement_type in ['unanimous', 'mostly_agree', 'ambiguous']:
        count = top5_agreement[agreement_type]
        pct = 100 * count / total_reads if total_reads else 0
        if agreement_type == 'unanimous':
            print(f"All top-5 hits same taxonomy: {count:,} reads ({pct:.2f}%)")
        elif agreement_type == 'mostly_agree':
            print(f"Top-5 hits split (2 taxa): {count:,} reads ({pct:.2f}%)")
        else:
            print(f"Ambiguous (3+ taxa in top-5): {count:,} reads ({pct:.2f}%)")
    print()

    # COVID-specific stats
    if covid_hits:
        all_covid_count = sum(1 for h in covid_hits if h['top5_all_covid'])
        print(f"COVID reads where all top-5 hits are SARS-CoV-2: {all_covid_count:,} / {len(covid_hits):,} "
              f"({100*all_covid_count/len(covid_hits):.1f}%)")
        print()

    # Show examples of ambiguous reads
    if ambiguous_reads:
        print("="*70)
        print(f"AMBIGUOUS READS (showing first 10 of {len(ambiguous_reads):,})")
        print("="*70)
        for qseqid, taxonomies, scinames in ambiguous_reads[:10]:
            print(f"Read: {qseqid}")
            print(f"  Top-5 taxonomies: {taxonomies}")
            print(f"  Top-5 organisms: {[s[:50] for s in scinames]}")
            print()
    print()

    # Write COVID TSV
    covid_tsv = f"{output_prefix}_covid_gene_mapping.tsv"
    with open(covid_tsv, 'w') as f:
        f.write("read_name\tmapped_length\tstart\tstop\tgene_name\tpercent_identity\tevalue\ttop5_all_covid\tnum_blast_hits\n")
        for hit in covid_hits:
            f.write(f"{hit['read_name']}\t{hit['mapped_length']}\t"
                   f"{hit['start']}\t{hit['end']}\t{hit['genes']}\t"
                   f"{hit['pident']:.2f}\t{hit['evalue']:.2e}\t"
                   f"{hit['top5_all_covid']}\t{hit['num_hits']}\n")

    print(f"Wrote {len(covid_hits):,} COVID reads to: {covid_tsv}")

    # Write non-COVID FASTAs by category
    for taxonomy, hits in non_covid_hits.items():
        fasta_out = f"{output_prefix}_{taxonomy}_reads.fasta"
        with open(fasta_out, 'w') as f:
            for read_id, sciname, pident, evalue in hits:
                if read_id in sequences:
                    # Clean up scientific name for header
                    sciname_clean = sciname.replace(' ', '_').replace('/', '_')
                    f.write(f">{read_id}|{sciname_clean}|pident={pident:.1f}|eval={evalue:.2e}\n")
                    f.write(f"{sequences[read_id]}\n")
        print(f"Wrote {len(hits):,} {taxonomy.upper()} reads to: {fasta_out}")

    print("\n" + "="*70)
    print("COMPLETE!")
    print("="*70)

if __name__ == '__main__':
    main()

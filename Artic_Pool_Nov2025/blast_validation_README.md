# BLAST Validation Workflow

Validates that mapped reads are truly SARS-CoV-2 and not human/bacterial contaminants.

## Quick Start

### Option 1: Full workflow (automated)
```bash
bash blast_validation_workflow.sh Handley_A2690_54317_Wilen_9_mappings.final.sam Wilen_9_blast
```

**WARNING**: Remote BLAST against nt is VERY slow (hours to days for thousands of reads).

### Option 2: Test with subset first
```bash
# Extract first 100 reads for testing
python3 extract_reads_for_blast.py Handley_A2690_54317_Wilen_9_mappings.final.sam \
    test_100_reads.fasta \
    --min-mapq 30 \
    --metadata-tsv test_100_metadata.tsv

# Take only first 100 sequences
head -200 test_100_reads.fasta > test_subset.fasta

# BLAST subset
blastn -query test_subset.fasta \
    -db nt \
    -remote \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \
    -max_target_seqs 1 \
    -out test_blast_results.txt

# Parse results
python3 parse_blast_results.py test_blast_results.txt test_subset.fasta test_100_metadata.tsv test_output
```

### Option 3: Step-by-step manual control

#### Step 1: Extract reads from SAM
```bash
python3 extract_reads_for_blast.py INPUT.sam OUTPUT.fasta \
    --min-mapq 30 \
    --metadata-tsv metadata.tsv
```

#### Step 2: BLAST against nt

**Remote BLAST** (uses NCBI servers, slow but no local database needed):
```bash
blastn -query reads.fasta \
    -db nt \
    -remote \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \
    -max_target_seqs 1 \
    -out blast_results.txt
```

**Local BLAST** (requires local nt database, much faster):
```bash
blastn -query reads.fasta \
    -db /path/to/nt \
    -num_threads 8 \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames' \
    -max_target_seqs 1 \
    -out blast_results.txt
```

#### Step 3: Parse BLAST results
```bash
python3 parse_blast_results.py blast_results.txt reads.fasta metadata.tsv output_prefix
```

## Output Files

### For COVID hits:
- `{prefix}_covid_gene_mapping.tsv` - TSV with columns:
  - read_name
  - mapped_length
  - start
  - stop
  - gene_name (e.g., "S", "N", "orf1ab", "5'_UTR")
  - percent_identity
  - evalue

### For non-COVID hits:
- `{prefix}_human_reads.fasta` - Reads matching human sequences
- `{prefix}_bacterial_reads.fasta` - Reads matching bacteria
- `{prefix}_viral_reads.fasta` - Reads matching other viruses
- `{prefix}_fungal_reads.fasta` - Reads matching fungi
- `{prefix}_other_reads.fasta` - Reads matching other organisms

Each FASTA header includes organism name, percent identity, and e-value.

## SARS-CoV-2 Genes Annotated

The script recognizes these genes (MN908947.3 coordinates):
- orf1ab (266-21555)
- S / Spike (21563-25384)
- ORF3a (25393-26220)
- E / Envelope (26245-26472)
- M / Membrane (26523-27191)
- ORF6 (27202-27387)
- ORF7a (27394-27759)
- ORF7b (27756-27887)
- ORF8 (27894-28259)
- N / Nucleocapsid (28274-29533)
- ORF10 (29558-29674)
- 5'_UTR (1-265)
- 3'_UTR (29534-29903)

## Performance Tips

1. **For large datasets**: Consider splitting FASTA into chunks and running parallel BLAST jobs
2. **For remote BLAST**: Use `-max_target_seqs 1` to get only top hit (faster)
3. **For testing**: Always test with a small subset first (100-1000 reads)
4. **Local nt database**: If available, much faster than remote BLAST

## Troubleshooting

**Remote BLAST timing out:**
- Split into smaller batches
- Use local BLAST if database is available
- Contact NCBI to ensure your IP isn't rate-limited

**No taxonomy information in BLAST results:**
- Ensure you're using the full format string with `staxids sscinames`
- Check that BLAST can access taxonomy database

**Reads classified incorrectly:**
- Check percent identity threshold (default: all hits accepted)
- Verify taxonomy classification logic in `parse_blast_results.py`
- Some sequences may legitimately match multiple organisms

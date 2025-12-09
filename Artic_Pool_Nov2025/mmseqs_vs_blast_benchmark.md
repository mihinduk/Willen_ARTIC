# MMseqs2 vs BLAST Benchmark Guide

Compare MMseqs2 and BLAST performance for validating SARS-CoV-2 reads.

## Quick Start

### Install MMseqs2 (if needed)
```bash
# Option 1: Conda
conda install -c bioconda mmseqs2

# Option 2: Module (if available on HTCF)
module load mmseqs2

# Check installation
mmseqs version
```

## Running Both Methods for Comparison

### 1. Submit Both Jobs

```bash
# Submit BLAST job
sbatch submit_blast_validation.sh

# Submit MMseqs2 job
sbatch submit_mmseqs_validation.sh

# Monitor both
squeue -u mihindu
```

### 2. Compare Performance

Track these metrics:
- **Runtime**: Check log files for total time
- **Memory usage**: Check SLURM accounting
- **Results**: Compare classification agreement

```bash
# After jobs complete:

# Check runtime
grep "Finished:" blast_validation_*.log mmseqs_validation_*.log

# Check memory
sacct -j <BLAST_JOB_ID> --format=JobID,Elapsed,MaxRSS
sacct -j <MMSEQS_JOB_ID> --format=JobID,Elapsed,MaxRSS

# Compare results
wc -l Wilen_9_blast_covid_gene_mapping.tsv Wilen_9_mmseqs_covid_gene_mapping.tsv
```

### 3. Compare Classification Results

```bash
# Count reads by category
echo "BLAST Results:"
grep -c "^>" Wilen_9_blast_covid_reads.fasta Wilen_9_blast_human_reads.fasta Wilen_9_blast_bacterial_reads.fasta 2>/dev/null

echo "MMseqs2 Results:"
grep -c "^>" Wilen_9_mmseqs_covid_reads.fasta Wilen_9_mmseqs_human_reads.fasta Wilen_9_mmseqs_bacterial_reads.fasta 2>/dev/null

# Check agreement for COVID reads
cut -f1 Wilen_9_blast_covid_gene_mapping.tsv | sort > blast_covid_reads.txt
cut -f1 Wilen_9_mmseqs_covid_gene_mapping.tsv | sort > mmseqs_covid_reads.txt

echo "Reads called COVID by both:"
comm -12 blast_covid_reads.txt mmseqs_covid_reads.txt | wc -l

echo "Only BLAST:"
comm -23 blast_covid_reads.txt mmseqs_covid_reads.txt | wc -l

echo "Only MMseqs2:"
comm -13 blast_covid_reads.txt mmseqs_covid_reads.txt | wc -l
```

## Expected Performance Differences

### BLAST
**Pros:**
- ✅ Standard, well-established tool
- ✅ Mature taxonomy integration
- ✅ Widely accepted for publications

**Cons:**
- ❌ Slower (6-24 hours for 700K reads)
- ❌ Higher memory usage

### MMseqs2
**Pros:**
- ✅ Much faster (30min - 2 hours for 700K reads)
- ✅ Lower memory footprint
- ✅ Comparable sensitivity to BLAST

**Cons:**
- ❌ Taxonomy annotation less mature
- ❌ Results may differ slightly from BLAST

## Typical Speedup

For this dataset (~700K reads):
- **BLAST**: 8-12 hours with 16 threads
- **MMseqs2**: 30min - 2 hours with 16 threads
- **Speedup**: ~10-20x faster

## Interpreting Differences

Small differences in classification are expected:

1. **Algorithm differences**: MMseqs2 uses different heuristics
2. **Sensitivity settings**: Default sensitivities differ slightly
3. **Taxonomy mapping**: Taxonomy assignment methods differ

**Action items if results differ significantly (>5%):**
- Adjust MMseqs2 sensitivity: `-s 7.5` (more sensitive) to `-s 4.0` (most sensitive)
- Check if one method has more "unknown" classifications
- Compare e-value distributions

## MMseqs2 Database Notes

The workflow will:
1. **Try to find existing MMseqs2 database** at `/ref/sahlab/data/nt/mmseqs_nt_db/nt`
2. **If not found, create one** from `/ref/sahlab/data/nt/nt.gz` (takes ~1-2 hours)
3. **Reuse the database** for subsequent runs

**Creating MMseqs2 nt database manually:**
```bash
# One-time setup (run as interactive job)
srun --mem=128G --cpus-per-task=16 --time=4:00:00 --pty bash

mkdir -p /ref/sahlab/data/nt/mmseqs_nt_db
gunzip -c /ref/sahlab/data/nt/nt.gz | mmseqs createdb stdin /ref/sahlab/data/nt/mmseqs_nt_db/nt --dbtype 2

# This creates ~500GB of index files and takes 1-2 hours
```

## Troubleshooting

**MMseqs2 not found:**
```bash
conda install -c bioconda mmseqs2
# OR
module load mmseqs2
```

**Out of memory:**
- Increase `--mem` in submit script
- Reduce `--threads`
- Split input FASTA into chunks

**No taxonomy information:**
- Results will still work but organisms will be "Unknown"
- Taxonomy is inferred from accession patterns
- SARS-CoV-2 detection is robust to this

## Recommendation

For this validation task:
1. **Run MMseqs2 first** (faster, good for initial results)
2. **Run BLAST for publication** (if reviewers prefer established tools)
3. **Report both** if results agree (shows robustness)

The two methods should give >95% agreement for a clean dataset like this.

# Processing Notes: MiSeq_i100_MD_327 (November 2025)

## Data Information

- **Run ID**: MiSeq_i100_MD_327
- **Source**: `/lts/sahlab/data4/DATA_DOWNLOADS_3/MiSeq_i100_MD_327`
- **Collaborator**: Craig Willen
- **Local Working Directory**: `/mnt/pathogen2/kathie/Artic_Pool_Nov2025`
- **Dropbox Location**: `Handley Lab Dropbox/virome/Willen_artic/Artic_Pool_Nov2025`

## Data Download via Globus

1. Follow the link in the notification email and log in to Globus
   - **Important**: Login as "Washington University in St. Louis" (NOT with ORCID)

2. In the search bar, enter: `HTCF@WUSTL`

3. In the path bar below HTCF@WUSTL, enter:
   ```
   /lts/sahlab/data4/DATA_DOWNLOADS_3/MiSeq_i100_MD_327
   ```

4. Select all files in the collection panel, click "Transfer or Sync", then "Start"

## Local Data Transfer

```bash
# Create directories
mkdir -p /mnt/pathogen2/kathie/Artic_Pool_Nov2025/raw
mkdir -p /mnt/pathogen2/kathie/Artic_Pool_Nov2025/artic_results

# Fetch FASTQs from HTCF
cd /mnt/pathogen2/kathie/Artic_Pool_Nov2025/raw
rsync -avh mihindu@login.htcf.wustl.edu:/lts/sahlab/data4/DATA_DOWNLOADS_3/MiSeq_i100_MD_327/*R1*gz .
rsync -avh mihindu@login.htcf.wustl.edu:/lts/sahlab/data4/DATA_DOWNLOADS_3/MiSeq_i100_MD_327/*R2*gz .
```

## Run Results

### Pangolin QC Summary
- **Failed samples**: 9 (based on `grep -c "fail"` on lineage output)
- **Status**: Processing stopped after Pangolin due to poor results

### Output Files Generated
- `2025_11_covid_combined.fa` - Combined consensus sequences
- `2025_11_covid.pangolin_lineage.csv` - Lineage assignments
- `2025_11_covid_combined_1000_N_or_less.fasta` - Filtered sequences (≤1000 Ns)

## Notes

### Quality Issues
The November 2025 run had significant quality issues:
- 9 samples failed Pangolin QC
- Pipeline stopped after step 9A (Pangolin) due to poor results
- Nextclade and downstream steps were not completed for this batch

### Follow-up
- Email sent to Jessica regarding security on 08/19/2025

## Resource Locations

| Resource | Path |
|----------|------|
| Reference FASTA | `/mnt/pathogen2/kathie/artic_resources/nCoV-2019.reference.fasta` |
| ARTIC Primer BED | `/mnt/pathogen2/kathie/artic_resources/ARTIC-V1.bed` |
| snpEff JAR | `/home/kathiem/snpEff/snpEff.jar` |
| SnpSift JAR | `/mnt/pathogen2/kathie/software/miniconda3/share/snpsift-4.3.1t-1/SnpSift.jar` |
| N-count script | `/Users/handley_lab/Resources/perl_scripts/N_count_and_clean.pl` |
| R analysis script | `2024_10_22_2024_covid.Rmd` |

## Pipeline Parameters Used

| Parameter | Value |
|-----------|-------|
| BWA threads | 12 |
| Samtools sort threads | 12 |
| LoFreq threads | 12 |
| LoFreq filter -v | 75 |
| mpileup depth (-d) | 1000 |
| N threshold for filtering | 1000 |

---
*Ars longa, vita brevis - Art is long, life is short - Hippocrates*

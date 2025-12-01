# ARTIC Pool November 2025 - SARS-CoV-2 Coverage Analysis

MiSeq i100 (2x250bp) sequencing analysis of SARS-CoV-2 samples from Wilen lab.

## Contents

### Scripts
- **`plot_coverage_map.py`** - Generate coverage maps from SAM files with gene annotations
  - Uses CIGAR strings to determine actual mapped positions (no artificial inflation)
  - Plots log10(depth+1) by default for better visualization
  - Outputs PNG and PDF formats
  - Shows 5' and 3' UTR regions

- **`calculate_median_mapped_length.py`** - Calculate median mapped read length after soft-clipping
  - Parses CIGAR strings to get true alignment lengths
  - Reports median, mean, min, max mapped lengths

### Reference Files
- **`nCoV-2019.reference.fasta`** - SARS-CoV-2 reference genome (MN908947.3)
- **`2025_11_covid_combined.fa`** - Combined sequences

### Metadata
- **`SARS_CoV2_plate_layout_submission_2025_11_04.xlsx`** - Sample plate layout
- **`2025_11_24_Download_of_ARTIC_MiSeq_i100_MD_327_processing.txt`** - Processing notes

### Coverage Maps
PNG and PDF files showing coverage across entire SARS-CoV-2 genome for each sample:
- Handley_A2682_54325_Wilen_1_5405_RECTUM
- Handley_A2683_54324_Wilen_2_5374
- Handley_A2684_54323_Wilen_3_5897_RECTUM
- Handley_A2685_54322_Wilen_4_COVID_5323_RECTUM
- Handley_A2686_54321_Wilen_5_LIINC_5104_RECTUM
- Handley_A2687_54320_Wilen_6_LIINC_5104_ILEUM
- Handley_A2688_54319_Wilen_7_FS12
- Handley_A2689_54318_Wilen_8_FS12
- Handley_A2690_54317_Wilen_9

### Presentation
- **`Willen_ARTIC_MiSeq_i100_MD_327_troubleshooting.key`** - Analysis and troubleshooting presentation

## Key Findings

**Critical Discovery**: Initial analysis using assumed 250bp read length artificially inflated coverage by 8-9x. Many reads had heavy soft-clipping (e.g., 221S30M CIGAR = only 30bp mapped out of 251bp total).

**Solution**: Updated analysis to parse CIGAR strings and use actual mapped positions, revealing true coverage depth.

## Usage

### Generate Coverage Maps
```bash
# From full SAM file (must have CIGAR strings)
python3 plot_coverage_map.py sample_mappings.nonsupplemental.sam --accession MN908947.3

# Adjust MAPQ threshold (default 30)
python3 plot_coverage_map.py sample.sam --accession MN908947.3 --min-mapq 60

# Use linear scale instead of log
python3 plot_coverage_map.py sample.sam --accession MN908947.3 --linear
```

### Calculate Median Mapped Length
```bash
# All reads
python3 calculate_median_mapped_length.py < sample.sam

# Filter by MAPQ >= 30
awk '$5 >= 30' sample.sam | python3 calculate_median_mapped_length.py
```

## Requirements

- Python 3.7+
- Biopython
- matplotlib
- numpy

```bash
pip install biopython matplotlib numpy
```

## Date
November-December 2025

---
*Veritas numquam perit - Truth never perishes*

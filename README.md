# Willen ARTIC SARS-CoV-2 Analysis Pipeline

This repository contains the ARTIC analysis workflow for Craig Willen's SARS-CoV-2 sequencing project.

## Overview

This pipeline processes Illumina MiSeq paired-end reads using the ARTIC protocol for SARS-CoV-2 genome assembly, variant calling, and lineage assignment.

## Pipeline Steps

1. **Data Download** - Transfer raw FASTQ files from HTCF via Globus
2. **Read Mapping** - Align reads to reference (MN908947.3) using BWA MEM
3. **Primer Trimming** - Remove ARTIC primer sequences using ivar trim
4. **Consensus Generation** - Generate consensus FASTA using ivar consensus
5. **Variant Calling** - Call variants with LoFreq (includes realignment and indel quality insertion)
6. **Variant Annotation** - Annotate variants using snpEff
7. **Lineage Assignment** - Assign Pangolin lineages and run Nextclade QC
8. **Quality Filtering** - Remove sequences with excessive Ns

## Reference Genome

- **Primary Reference**: MN908947.3 (ARTIC default, Wuhan-Hu-1)
- **snpEff Reference**: NC_045512.2 (for CZI pipeline compatibility)
- **Primer Scheme**: ARTIC V1

## Directory Structure

```
project_root/
├── raw/                    # Raw FASTQ files (R1 and R2)
├── artic_results/          # All output files
│   ├── *.sorted.bam        # Initial aligned BAMs
│   ├── *.primertrim.sorted.bam  # Primer-trimmed BAMs
│   ├── *.consensus.fa      # Individual consensus sequences
│   ├── *.lofreq.final.bam  # LoFreq-processed BAMs
│   ├── *_vars.filt.vcf     # Filtered variant calls
│   ├── *.snpEFF.ann.vcf    # Annotated variants
│   └── *.snpsift.txt       # Extracted variant fields
└── scripts/                # Processing scripts
```

## Dependencies

### Software
- BWA (alignment)
- samtools (BAM manipulation)
- ivar (primer trimming, consensus calling)
- LoFreq (variant calling)
- snpEff (variant annotation)
- SnpSift (variant field extraction)
- Pangolin (lineage assignment)

### Conda Installation
```bash
conda install -c bioconda bwa samtools ivar lofreq pangolin
conda install -c conda-forge ncurses
conda install -c "bioconda/label/cf201901" snpsift
```

## Usage

### Quick Start
```bash
# Set up directories
mkdir -p /path/to/project/raw
mkdir -p /path/to/project/artic_results

# Set variables
IN=./raw
OUT=./artic_results
REF=$OUT/nCoV-2019.reference.fasta

# Run the pipeline
bash scripts/artic_processing.sh
```

### Running Individual Steps

See `scripts/artic_processing.sh` for the complete pipeline, or run individual steps as documented in `docs/processing_notes.md`.

## Outputs

### Key Output Files
- `*_combined.fa` - Combined consensus sequences for all samples
- `*.pangolin_lineage.csv` - Pangolin lineage assignments
- `*.snpEFF.ann.tsv` - Annotated variants in tabular format
- `*.snpsift.txt` - Extracted variant fields (POS, REF, ALT, DP4, EFFECT, GENE, HGVS_P)

### Quality Metrics
- Pangolin QC status indicates pass/fail for each sample
- N-count filtering removes sequences with >1000 Ns (configurable)

## Post-Processing

### Nextclade Analysis
Upload combined consensus FASTA to [Nextclade](https://clades.nextstrain.org/) for additional QC and clade assignment.

### Data Visualization
R script `2024_10_22_2024_covid.Rmd` joins lineage data and produces final filtered output.

## Collaborators

- **Lab**: Handley Lab (sahlab), Washington University School of Medicine
- **Collaborator**: Craig Willen
- **Data Source**: MiSeq_i100_MD_327

## Processing Logs

See `docs/` directory for detailed processing notes and run-specific information.

## References

- [ARTIC Network](https://artic.network/)
- [Pangolin](https://cov-lineages.org/pangolin.html)
- [Nextclade](https://clades.nextstrain.org/)
- [ivar](https://github.com/andersen-lab/ivar)
- [LoFreq](https://csb5.github.io/lofreq/)

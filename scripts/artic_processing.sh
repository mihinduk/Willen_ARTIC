#!/bin/bash
#===============================================================================
# ARTIC SARS-CoV-2 Processing Pipeline
# 
# Description: Process Illumina MiSeq reads using ARTIC protocol
# Author: Handley Lab (sahlab)
# Date: November 2025
#===============================================================================

set -e  # Exit on error

#-------------------------------------------------------------------------------
# Configuration - MODIFY THESE FOR EACH RUN
#-------------------------------------------------------------------------------
PROJECT_DIR="/mnt/pathogen2/kathie/Artic_Pool_Nov2025"
RUN_NAME="2025_11_covid"
THREADS=12
N_THRESHOLD=1000  # Max Ns allowed in consensus sequences

# Derived paths
IN="${PROJECT_DIR}/raw"
OUT="${PROJECT_DIR}/artic_results"
REF="${OUT}/nCoV-2019.reference.fasta"

# Resource paths
ARTIC_RESOURCES="/mnt/pathogen2/kathie/artic_resources"
SNPEFF_JAR="/home/kathiem/snpEff/snpEff.jar"
SNPSIFT_JAR="/mnt/pathogen2/kathie/software/miniconda3/share/snpsift-4.3.1t-1/SnpSift.jar"
N_COUNT_SCRIPT="${OUT}/N_count_and_clean.pl"

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------
log_step() {
    echo ""
    echo "=================================================================="
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
    echo "=================================================================="
    echo ""
}

check_file() {
    if [[ ! -f "$1" ]]; then
        echo "ERROR: Required file not found: $1"
        exit 1
    fi
}

#-------------------------------------------------------------------------------
# Step 0: Setup directories
#-------------------------------------------------------------------------------
setup_directories() {
    log_step "Setting up directories"
    mkdir -p "$IN"
    mkdir -p "$OUT"
}

#-------------------------------------------------------------------------------
# Step 1: Fetch reference and primer files
#-------------------------------------------------------------------------------
fetch_resources() {
    log_step "Fetching reference and primer files"
    cd "$OUT"
    rsync -avh "${ARTIC_RESOURCES}/nCoV-2019.reference.fasta"* .
    rsync -avh "${ARTIC_RESOURCES}/ARTIC-V1.bed" .
    cd "$PROJECT_DIR"
}

#-------------------------------------------------------------------------------
# Step 2: Map reads to reference
#-------------------------------------------------------------------------------
map_reads() {
    log_step "Mapping reads to reference with BWA MEM"
    
    for i in "$IN"/*_R1_001.fastq.gz; do
        sample=$(basename "$i" _R1_001.fastq.gz)
        echo "Mapping: $sample"
        
        bwa mem -t "$THREADS" "$REF" \
            "${IN}/${sample}_R1_001.fastq.gz" \
            "${IN}/${sample}_R2_001.fastq.gz" \
            | samtools sort \
            | samtools view -F 4 -o "${OUT}/${sample}.sorted.bam"
    done
}

#-------------------------------------------------------------------------------
# Step 3: Trim primers with ivar
#-------------------------------------------------------------------------------
trim_primers() {
    log_step "Trimming primers with ivar"
    
    for i in "$OUT"/*.sorted.bam; do
        sample=$(basename "$i" .sorted.bam)
        echo "Primer trim: $sample"
        
        ivar trim -e \
            -i "${OUT}/${sample}.sorted.bam" \
            -b "${OUT}/ARTIC-V1.bed" \
            -p "${OUT}/${sample}.primertrim"
    done
}

#-------------------------------------------------------------------------------
# Step 4: Re-sort primer-trimmed BAMs
#-------------------------------------------------------------------------------
resort_bams() {
    log_step "Re-sorting primer-trimmed BAMs"
    
    for i in "$OUT"/*.primertrim.bam; do
        sample=$(basename "$i" .primertrim.bam)
        echo "Resort: $sample"
        
        samtools sort "${OUT}/${sample}.primertrim.bam" \
            -o "${OUT}/${sample}.primertrim.sorted.bam"
    done
}

#-------------------------------------------------------------------------------
# Step 5: Generate consensus sequences
#-------------------------------------------------------------------------------
generate_consensus() {
    log_step "Generating consensus sequences with ivar"
    
    for i in "$OUT"/*.primertrim.sorted.bam; do
        sample=$(basename "$i" .primertrim.sorted.bam)
        echo "Consensus: $sample"
        
        samtools mpileup -A -d 1000 -B -Q 0 \
            --reference "$REF" \
            "${OUT}/${sample}.primertrim.sorted.bam" \
            | ivar consensus -p "${OUT}/${sample}.consensus" -n N
    done
}

#-------------------------------------------------------------------------------
# Step 6: LoFreq variant calling pipeline
#-------------------------------------------------------------------------------
lofreq_realign() {
    log_step "LoFreq: Viterbi realignment"
    
    for i in "$OUT"/*.primertrim.sorted.bam; do
        sample=$(basename "$i" .primertrim.sorted.bam)
        echo "LoFreq realign: $sample"
        
        lofreq viterbi -f "$REF" "${OUT}/${sample}.primertrim.sorted.bam" \
            | samtools sort -n --threads "$THREADS" \
            -o "${OUT}/${sample}.lofreq.realign.bam"
    done
}

lofreq_indelqual() {
    log_step "LoFreq: Insert indel qualities"
    
    for i in "$OUT"/*.lofreq.realign.bam; do
        sample=$(basename "$i" .lofreq.realign.bam)
        echo "LoFreq indelqual step 1: $sample"
        
        lofreq indelqual --dindel -f "$REF" "${OUT}/${sample}.lofreq.realign.bam" \
            | samtools sort --threads "$THREADS" \
            -o "${OUT}/${sample}.lofreq.indel.bam"
    done
    
    for i in "$OUT"/*.lofreq.realign.bam; do
        sample=$(basename "$i" .lofreq.realign.bam)
        echo "LoFreq alnqual step 2: $sample"
        
        lofreq alnqual -b "${OUT}/${sample}.lofreq.indel.bam" "$REF" \
            > "${OUT}/${sample}.lofreq.final.bam"
    done
    
    for i in "$OUT"/*.lofreq.final.bam; do
        sample=$(basename "$i" .lofreq.final.bam)
        echo "Indexing: $sample"
        
        samtools index "${OUT}/${sample}.lofreq.final.bam"
    done
}

lofreq_call_variants() {
    log_step "LoFreq: Call variants"
    
    for i in "$OUT"/*.lofreq.final.bam; do
        sample=$(basename "$i" .lofreq.final.bam)
        echo "Call variants: $sample"
        
        lofreq call-parallel --pp-threads "$THREADS" \
            --force-overwrite \
            --no-default-filter \
            --call-indels \
            -f "$REF" \
            -o "${OUT}/${sample}_vars.vcf" \
            "${OUT}/${sample}.lofreq.final.bam"
    done
}

lofreq_filter() {
    log_step "LoFreq: Filter variants"
    
    for i in "$OUT"/*_vars.vcf; do
        sample=$(basename "$i" _vars.vcf)
        echo "Filtering: $sample"
        
        lofreq filter -i "$i" \
            -o "${OUT}/${sample}_vars.filt.vcf" \
            -v 75
    done
}

#-------------------------------------------------------------------------------
# Step 7: snpEff annotation
#-------------------------------------------------------------------------------
run_snpeff() {
    log_step "Running snpEff annotation"
    
    for i in "$OUT"/*_vars.filt.vcf; do
        sample=$(basename "$i" _vars.filt.vcf)
        echo "snpEff: $sample"
        
        java -jar "$SNPEFF_JAR" MN908947.3 "$i" \
            -s "${OUT}/${sample}_summary.html" \
            > "${OUT}/${sample}.snpEFF.ann.vcf"
        
        # Extract annotation fields to TSV
        grep -v "^##" "${OUT}/${sample}.snpEFF.ann.vcf" | \
            tail -n+2 | \
            cut -f8 | \
            sed 's/|/\t/g' | \
            cut -f1-16 | \
            sed '1i INFO\tEFFECT\tPUTATIVE_IMPACT\tGENE_NAME\tGENE_ID\tFEATURE_TYPE\tFEATURE_ID\tTRANSCRIPT_TYPE\tEXON_INTRON_RANK\tHGVSc\tHGVSp\tcDNA_POSITION_AND_LENGTH\tCDS_POSITION_AND_LENGTH\tPROTEIN_POSITION_AND_LENGTH\tDISTANCE_TO_FEATURE\tERROR' \
            > "${OUT}/${sample}.snpEFF.ann.tmp"
        
        grep -v "^##" "${OUT}/${sample}.snpEFF.ann.vcf" | \
            cut -f1-7 > "${OUT}/${sample}.ann.base.vcf"
        
        paste "${OUT}/${sample}.ann.base.vcf" "${OUT}/${sample}.snpEFF.ann.tmp" \
            > "${OUT}/${sample}.snpEFF.ann.tsv"
        
        rm "${OUT}/${sample}.snpEFF.ann.tmp"
        rm "${OUT}/${sample}.ann.base.vcf"
    done
}

#-------------------------------------------------------------------------------
# Step 8: Pangolin lineage assignment
#-------------------------------------------------------------------------------
run_pangolin() {
    log_step "Running Pangolin lineage assignment"
    
    # Combine consensus sequences
    cat "$OUT"/*consensus.fa > "${OUT}/${RUN_NAME}_combined.fa"
    
    # Run Pangolin
    pangolin "${OUT}/${RUN_NAME}_combined.fa" \
        --outfile "${OUT}/${RUN_NAME}.pangolin_lineage.csv"
    
    # Report failures
    fail_count=$(grep -c "fail" "${OUT}/${RUN_NAME}.pangolin_lineage.csv" || echo "0")
    echo "Samples failed Pangolin QC: $fail_count"
}

#-------------------------------------------------------------------------------
# Step 9: VCF conversion for CZI pipeline compatibility
#-------------------------------------------------------------------------------
convert_vcf_reference() {
    log_step "Converting VCF reference names (MN908947.3 -> NC_045512.2)"
    
    for i in "$OUT"/*_vars.filt.vcf; do
        sample=$(basename "$i" _vars.filt.vcf)
        echo "Converting: $sample"
        
        sed 's/MN908947.3/NC_045512.2/g' "${OUT}/${sample}_vars.filt.vcf" \
            > "${OUT}/${sample}.converted.vcf"
    done
}

run_snpeff_nc() {
    log_step "Running snpEff with NC_045512.2 reference"
    
    for i in "$OUT"/*.converted.vcf; do
        sample=$(basename "$i" .converted.vcf)
        echo "snpEff NC: $sample"
        
        java -Xmx8g -jar "$SNPEFF_JAR" NC_045512.2 \
            "${OUT}/${sample}.converted.vcf" \
            > "${OUT}/${sample}.converted.ann.vcf"
        
        mv snpEff_genes.txt "${OUT}/${sample}_snpEff_genes.txt"
        mv snpEff_summary.html "${OUT}/${sample}_snpEff_summary.html"
    done
}

run_snpsift() {
    log_step "Running SnpSift field extraction"
    
    for i in "$OUT"/*.converted.ann.vcf; do
        sample=$(basename "$i" .converted.ann.vcf)
        echo "SnpSift: $sample"
        
        java -jar "$SNPSIFT_JAR" extractFields \
            "${OUT}/${sample}.converted.ann.vcf" \
            POS REF ALT DP4 \
            "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].HGVS_P" \
            > "${OUT}/${sample}.snpsift.txt"
    done
}

#-------------------------------------------------------------------------------
# Step 10: Quality filtering
#-------------------------------------------------------------------------------
filter_by_n_count() {
    log_step "Filtering sequences by N count (threshold: $N_THRESHOLD)"
    
    check_file "$N_COUNT_SCRIPT"
    
    cd "$OUT"
    perl "$N_COUNT_SCRIPT" "${RUN_NAME}_combined.fa" "$N_THRESHOLD"
    cd "$PROJECT_DIR"
}

#-------------------------------------------------------------------------------
# Main Pipeline Execution
#-------------------------------------------------------------------------------
run_full_pipeline() {
    log_step "Starting ARTIC Processing Pipeline"
    echo "Project: $PROJECT_DIR"
    echo "Run name: $RUN_NAME"
    echo "Threads: $THREADS"
    
    setup_directories
    fetch_resources
    map_reads
    trim_primers
    resort_bams
    generate_consensus
    lofreq_realign
    lofreq_indelqual
    lofreq_call_variants
    lofreq_filter
    run_snpeff
    run_pangolin
    convert_vcf_reference
    run_snpeff_nc
    run_snpsift
    filter_by_n_count
    
    log_step "Pipeline Complete!"
    echo "Results in: $OUT"
}

#-------------------------------------------------------------------------------
# Command Line Interface
#-------------------------------------------------------------------------------
usage() {
    echo "Usage: $0 [OPTION]"
    echo ""
    echo "Options:"
    echo "  --full          Run full pipeline"
    echo "  --map           Run mapping only"
    echo "  --trim          Run primer trimming only"
    echo "  --consensus     Generate consensus only"
    echo "  --variants      Run variant calling only"
    echo "  --annotate      Run annotation only"
    echo "  --pangolin      Run Pangolin only"
    echo "  --help          Show this help message"
    echo ""
}

case "${1:-}" in
    --full)
        run_full_pipeline
        ;;
    --map)
        map_reads
        ;;
    --trim)
        trim_primers
        resort_bams
        ;;
    --consensus)
        generate_consensus
        ;;
    --variants)
        lofreq_realign
        lofreq_indelqual
        lofreq_call_variants
        lofreq_filter
        ;;
    --annotate)
        run_snpeff
        convert_vcf_reference
        run_snpeff_nc
        run_snpsift
        ;;
    --pangolin)
        run_pangolin
        ;;
    --help|"")
        usage
        ;;
    *)
        echo "Unknown option: $1"
        usage
        exit 1
        ;;
esac

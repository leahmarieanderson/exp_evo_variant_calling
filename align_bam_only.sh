#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Leah/labmeeting_250613
#$ -o /net/dunham/vol2/Leah/labmeeting_250613/outputs/
#$ -e /net/dunham/vol2/Leah/labmeeting_250613/errors/
#$ -l mfree=8G
#$ -l h_rt=36:0:0


## Alignment-only pipeline for YEvo data
## Modified from original align.sh script
## Performs FastQC, BWA alignment, duplicate marking, and creates final BAM file

module load modules modules-init modules-gs
module load zlib/1.3.1
module load bwa/0.7.17
module load htslib/1.19
module load samtools/1.19
module load picard/3.1.1
module load python/3.12.1 numpy biopython lofreq/2.1.5-18
module load java/1.17
module load GATK/4.5.0.0
module load perl/5.38.2
module load VCFtools/0.1.16-20
module load bcftools/1.20
module load bedtools/2.31.1
module load freebayes/1.3.6
module load fastqc/0.12.1


FOLDER=fastq
SAMPLE=$1 # Passed sample prefix (ex: Sample-01)
DIR=/net/dunham/vol2/Leah/labmeeting_250613
WORKDIR=${DIR}/WorkDirectory # Where files will be created
SEQDIR=${DIR}/${FOLDER} # Location of Fastqs
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
SEQID=leah-labmeeting # Project name and date for bam header
REF=${DIR}/exp_evo_variant_calling/genomes/sacCer3.fasta # Reference genome

# Sets up folder structure
mkdir -p ${WORKDIR}/${SAMPLE}
cd ${WORKDIR}/${SAMPLE}

# Perform FastQC checks on our samples
(>&2 echo ***FASTQC on SAMPLE ***)

# Check for different FASTQ file patterns
if ls ${SEQDIR}/${SAMPLE}_*R1*.fastq.gz 1> /dev/null 2>&1; then
    # Paired-end .fastq.gz files
    FASTQ_TYPE="paired_gz"
    fastqc ${SEQDIR}/${SAMPLE}_*R1*.fastq.gz -o ${WORKDIR}/${SAMPLE}/
    fastqc ${SEQDIR}/${SAMPLE}_*R2*.fastq.gz -o ${WORKDIR}/${SAMPLE}/
    # remove the zip file since the html file has everything we need to know
    rm ${SAMPLE}_*R1*fastqc.zip
    rm ${SAMPLE}_*R2*fastqc.zip
elif ls ${SEQDIR}/${SAMPLE}*.fastq.zip 1> /dev/null 2>&1; then
    # Single-end .fastq.zip files - need to extract first
    FASTQ_TYPE="single_zip"
    # Extract the zip file
    unzip -o ${SEQDIR}/${SAMPLE}*.fastq.zip -d ${SEQDIR}/
    # Find the extracted fastq file
    EXTRACTED_FASTQ=$(find ${SEQDIR}/ -name "${SAMPLE}*.fastq" -type f | head -1)
    if [ -z "$EXTRACTED_FASTQ" ]; then
        echo "Error: Could not find extracted FASTQ file for ${SAMPLE}"
        exit 1
    fi
    fastqc ${EXTRACTED_FASTQ} -o ${WORKDIR}/${SAMPLE}/
    # remove the fastqc zip file
    rm ${SAMPLE}*fastqc.zip
else
    echo "Error: Could not find FASTQ files for sample ${SAMPLE}"
    echo "Looking for either ${SAMPLE}_*R1*.fastq.gz and ${SAMPLE}_*R2*.fastq.gz"
    echo "or ${SAMPLE}*.fastq.zip"
    exit 1
fi

# Align reads with bwa
(>&2 echo ***BWA - mem -R***)

if [ "$FASTQ_TYPE" = "paired_gz" ]; then
    # Paired-end alignment
    bwa mem -R "@RG\tID:${SEQID}\tSM:${SAMPLE}\tLB:1" ${REF} ${SEQDIR}/${SAMPLE}_*R1*.fastq.gz ${SEQDIR}/${SAMPLE}_*R2*.fastq.gz > ${SAMPLE}_R1R2.sam
elif [ "$FASTQ_TYPE" = "single_zip" ]; then
    # Single-end alignment
    bwa mem -R "@RG\tID:${SEQID}\tSM:${SAMPLE}\tLB:1" ${REF} ${EXTRACTED_FASTQ} > ${SAMPLE}_R1R2.sam
fi

mkdir -p dup_metrics

(>&2 echo ***GATK4 - MarkDuplicatesSpark and Sort***)
# Previous align script had duplicates removed, this just has them marked down but not removed
gatk MarkDuplicatesSpark \
         -I ${SAMPLE}_R1R2.sam \
         -M dup_metrics/${SAMPLE}_dup_metrics.txt \
         -O ${SAMPLE}_R1R2_MD.sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${WORKDIR}/${SAMPLE}/${SAMPLE}_R1R2_MD.sort.bam

# Remove intermediate files
rm ${SAMPLE}_R1R2.sam

(>&2 echo ***ALIGNMENT COMPLETE - BAM file created: ${SAMPLE}_R1R2_MD.sort.bam***)
#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Leah/plasmidsaurus_wgs
#$ -o /net/dunham/vol2/Leah/plasmidsaurus_wgs/outputs/
#$ -e /net/dunham/vol2/Leah/plasmidsaurus_wgs/errors/
#$ -l mfree=8G
#$ -l h_rt=36:0:0


## SNP calling and alignment pipeline for YEvo data
## Chris Large and Caiti S. Heil. Modified for Bryce Taylor and Ryan Skophammer
## Uses the recommended SNP calling pipeline from Samtools
## Then filters based on the ancestral sequence

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


SAMPLE=fks1del # Passed sample prefix (ex: Sample-01)
DIR=/net/dunham/vol2/Leah/plasmidsaurus_wgs
WORKDIR=${DIR}
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
SEQID=fks1del # Project name and date for bam header
REF=${DIR}/exp_evo_variant_calling/genomes/sacCer3.fasta # Reference genome
ANNOTATE=${SCRIPTS}/genomes # Location of custom annotation scripts

# Set up working directory
cd ${WORKDIR}

# Extract ZIP file and find FASTQ
unzip -o "fks1del.fastq.zip"
EXTRACTED_FASTQ=$(find . -name "fks1del*.fastq" -type f | head -1)
if [ -z "$EXTRACTED_FASTQ" ]; then
    echo "Error: Could not find extracted FASTQ file for ${SAMPLE}"
    exit 1
fi

# Single-end alignment
bwa mem -R "@RG\tID:${SEQID}\tSM:${SAMPLE}\tLB:1" ${REF} ${EXTRACTED_FASTQ} > ${SAMPLE}_R1R2.sam


mkdir -p dup_metrics

(>&2 echo ***GATK4 - MarkDuplicatesSpark and Sort***)
# Previous align script had duplicates removed, this just has them marked down but not removed
gatk MarkDuplicatesSpark \
         -I ${SAMPLE}_R1R2.sam \
         -M dup_metrics/${SAMPLE}_dup_metrics.txt \
         -O ${SAMPLE}_R1R2_MD.sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${SAMPLE}_R1R2_MD.sort.bam

# Remove intermediate files
rm ${SAMPLE}_R1R2.sam


#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Zilong/updating_pipeline_2024
#$ -o /net/dunham/vol2/Zilong/updating_pipeline_2024/outputs/
#$ -e /net/dunham/vol2/Zilong/updating_pipeline_2024/errors/
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
module load GATK/3.7
module load python/3.12.1 numpy biopython
module load perl/5.38.2
module load VCFtools/0.1.16-20
module load bcftools/1.20
module load bedtools/2.31.1
module load freebayes/1.3.6


FOLDER=fastq
ANC=$1
DIR=/net/dunham/vol2/Zilong/updating_pipeline_2024
WORKDIR=${DIR}/WorkDirectory # Where files will be created
SEQDIR=${DIR}/${FOLDER} # Location of Fastqs
SEQID=leah_freeze_evolution # Project name and date for bam header
REF=${DIR}/genomes/sacCer3.fasta # Reference genome

# Sets up folder structure
mkdir -p ${WORKDIR}/${ANC}
cd ${WORKDIR}/${ANC}
#
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem -R '@RG\tID:'${SEQID}'\tSM:'${ANC}'\tLB:1' ${REF} ${SEQDIR}/${ANC}_*R1*.fastq.gz ${SEQDIR}/${ANC}_*R2*.fastq.gz > ${ANC}_R1R2.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${ANC}_R1R2.sam -o ${ANC}_R1R2.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${ANC}_R1R2.bam -o ${ANC}_R1R2_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${ANC}_R1R2_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${WORKDIR}/${ANC}/${ANC}_R1R2_sort.bam

# Remove intermediate files
rm ${ANC}_R1R2.sam
rm ${ANC}_R1R2.bam

cd ${WORKDIR}/${ANC}

mkdir -p dup_metrics

(>&2 echo ***Picard - MarkDuplicates***)
java -Xmx2g -jar $PICARD_DIR/picard.jar MarkDuplicates \
        INPUT=${ANC}_R1R2_sort.bam \
        OUTPUT=${ANC}_comb_R1R2.MD.bam \
        METRICS_FILE=dup_metrics/${ANC}_comb_R1R2.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
java -Xmx2g -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups \
        I=${ANC}_comb_R1R2.MD.bam \
        O=${ANC}_comb_R1R2.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${ANC} \
        VALIDATION_STRINGENCY=LENIENT

module unload picard/3.1.1
module load java/1.8.0

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${ANC}_comb_R1R2.RG.MD.bam \
        -o ${ANC}_comb_R1R2.RG.MD.sort.bam
samtools index ${ANC}_comb_R1R2.RG.MD.sort.bam

(>&2 echo ***GATK - RealingerTargetCreator***)
#GATK Realinger
java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R ${REF} \
        -I ${ANC}_comb_R1R2.RG.MD.sort.bam \
        -o ${ANC}_comb_R1R2.bam.intervals
java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R ${REF} \
        -I ${ANC}_comb_R1R2.RG.MD.sort.bam \
        -targetIntervals ${ANC}_comb_R1R2.bam.intervals \
        -o ${ANC}_comb_R1R2.RG.MD.realign.bam

samtools sort ${ANC}_comb_R1R2.RG.MD.realign.bam \
        -o ${ANC}_comb_R1R2.RG.MD.realign.sort.bam
samtools index ${ANC}_comb_R1R2.RG.MD.realign.sort.bam

cd ${WORKDIR}/${ANC}/

# Samtools pileup variant calling
(>&2 echo ***BCFtools - Pileup***)
bcftools mpileup --ignore-RG -Ou -ABf ${REF} ${ANC}_comb_R1R2.RG.MD.realign.sort.bam | bcftools call -vmO v -o ${ANC}_samtools_AB.vcf

# Freebayes with a lot of arguments for population calling
freebayes -f ${REF} \
        --pooled-continuous --report-genotype-likelihood-max --allele-balance-priors-off --min-alternate-fraction 0.1 \
        ${ANC}_comb_R1R2.RG.MD.realign.sort.bam > ${ANC}_freebayes_BCBio.vcf

# Remove intermediates
rm ${ANC}_comb_R1R2.bam.intervals
rm ${ANC}_comb_R1R2.MD.bam
rm ${ANC}_comb_R1R2.RG.MD.bam
rm ${ANC}_comb_R1R2.RG.MD.realign.bam
rm ${ANC}_comb_R1R2.RG.MD.realign.bai
rm ${ANC}_R1R2_sort.bam
rm ${ANC}_R1R2_sort.bam.bai
rm ${ANC}_comb_R1R2.RG.MD.sort.bam
rm ${ANC}_comb_R1R2.RG.MD.sort.bam.bai

mkdir bams

mv ${ANC}_comb_R1R2.RG.MD.realign.sort.bam bams
mv ${ANC}_comb_R1R2.RG.MD.realign.sort.bam.bai bams

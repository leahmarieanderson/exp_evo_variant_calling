#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Zilong/updating_pipeline_2024
#$ -o /net/dunham/vol2/Zilong/updating_pipeline_2024/outputs/
#$ -e /net/dunham/vol2/Zilong/updating_pipeline_2024/errors/
#$ -l mfree=8G
#$ -l h_rt=36:0:0

#set -e
#set -u

## SNP calling and alignment pipeline for YEvo data
## Chris Large and Caiti S. Heil. Modified for Bryce Taylor and Ryan Skophammer
## Uses the recommended SNP calling pipeline from Samtools
## Then filters based on the ancestral sequence

module load modules modules-init modules-gs
module load zlib/1.2.13_compat
module load bwa/0.7.15
module load htslib/1.18
module load samtools/1.14
module load picard/3.1.1
module load GATK/3.7
module load python/3.12.1 numpy biopython
module load perl/5.26.3
module load VCFtools/0.1.16-20
module load bcftools/1.19
module load bedtools/2.25.0
module load freebayes/1.3.6


FOLDER=fastq
SAMPLE=$1 # Passed sample prefix (ex: Sample-01)
ANC=$2
DIR=/net/dunham/vol2/Zilong/updating_pipeline_2024
WORKDIR=${DIR}/WorkDirectory # Where files will be created
SEQDIR=${DIR}/${FOLDER} # Location of Fastqs
SEQID=leah_freeze_evolution # Project name and date for bam header
REF=/net/dunham/vol2/Zilong/updating_pipeline_2024/genomes/sacCer3.fasta # Reference genome
ANNOTATE=/net/dunham/vol2/Cris_L/ReferenceGenome/S288C_reference_genome_R64-1-1_20110203 # Location of custom annotation scripts
SCRIPTS=/net/dunham/vol2/Zilong/updating_pipeline_2024/exp_evo_variant_calling/ # Path of annotation_final.py directory
ANCBAM=${WORKDIR}/${ANC}/${ANC}_comb_R1R2.RG.MD.realign.sort.bam
VCFDIR=${WORKDIR}/${ANC}/

# Sets up folder structure
mkdir -p ${WORKDIR}/${SAMPLE}
cd ${WORKDIR}/${SAMPLE}
#
# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem -R '@RG\tID:'${SEQID}'\tSM:'${SAMPLE}'\tLB:1' ${REF} ${SEQDIR}/${SAMPLE}_*R1*.fastq.gz ${SEQDIR}/${SAMPLE}_*R2*.fastq.gz > ${SAMPLE}_R1R2.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${SAMPLE}_R1R2.sam -o ${SAMPLE}_R1R2.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${SAMPLE}_R1R2.bam -o ${SAMPLE}_R1R2_sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${SAMPLE}_R1R2_sort.bam

# Print stats on how well the alignment worked
(>&2 echo ***Samtools - Flagstat***)
samtools flagstat ${WORKDIR}/${SAMPLE}/${SAMPLE}_R1R2_sort.bam

# Remove intermediate files
rm ${SAMPLE}_R1R2.sam
rm ${SAMPLE}_R1R2.bam

cd ${WORKDIR}/${SAMPLE}

mkdir -p dup_metrics

(>&2 echo ***Picard - MarkDuplicates***)
java -Xmx2g -jar $PICARD_DIR/picard.jar MarkDuplicates \
        INPUT=${SAMPLE}_R1R2_sort.bam \
        OUTPUT=${SAMPLE}_comb_R1R2.MD.bam \
        METRICS_FILE=dup_metrics/${SAMPLE}_comb_R1R2.sort_dup_metrics \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

(>&2 echo ***Picard - AddOrReplaceReadGroups***)
#Add or replace read groups needs to happen before GATK
java -Xmx2g -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups \
        I=${SAMPLE}_comb_R1R2.MD.bam \
        O=${SAMPLE}_comb_R1R2.RG.MD.bam \
        RGID=${SEQID} \
        RGLB=1 \
        RGPU=1 \
        RGPL=illumina \
        RGSM=${SAMPLE} \
        VALIDATION_STRINGENCY=LENIENT

module unload picard/3.1.1
module load java/1.8.0

(>&2 echo ***Samtools - Sort and Index***)
samtools sort ${SAMPLE}_comb_R1R2.RG.MD.bam \
        -o ${SAMPLE}_comb_R1R2.RG.MD.sort.bam
samtools index ${SAMPLE}_comb_R1R2.RG.MD.sort.bam

(>&2 echo ***GATK - RealingerTargetCreator***)
#GATK Realinger
java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R ${REF} \
        -I ${SAMPLE}_comb_R1R2.RG.MD.sort.bam \
        -o ${SAMPLE}_comb_R1R2.bam.intervals
java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R ${REF} \
        -I ${SAMPLE}_comb_R1R2.RG.MD.sort.bam \
        -targetIntervals ${SAMPLE}_comb_R1R2.bam.intervals \
        -o ${SAMPLE}_comb_R1R2.RG.MD.realign.bam

samtools sort ${SAMPLE}_comb_R1R2.RG.MD.realign.bam \
        -o ${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam
samtools index ${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam

cd ${WORKDIR}/${SAMPLE}/

# Samtools pileup variant calling
(>&2 echo ***BCFtools - Pileup***)
bcftools mpileup --ignore-RG -Ou -ABf ${REF} ${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam | bcftools call -vmO v -o ${SAMPLE}_samtools_AB.vcf

# Freebayes with a lot of arguments for population calling
freebayes -f ${REF} \
        --pooled-discrete --pooled-continuous --report-genotype-likelihood-max --allele-balance-priors-off --min-alternate-fraction 0.1 \
        ${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam > ${SAMPLE}_freebayes_BCBio.vcf

# Filters samtools by ancestor
(>2 echo ***Bedtools - Intersect***)
bedtools intersect -v -header \
        -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_samtools_AB.vcf \
        -b ${VCFDIR}/${ANC}_samtools_AB.vcf \
        > ${WORKDIR}/${SAMPLE}/${SAMPLE}_samtools_AB_AncFiltered.vcf

# Filters freebayes by ancestor 
bedtools intersect -v -header \
        -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio.vcf \
        -b ${VCFDIR}/${ANC}_freebayes_BCBio.vcf \
        > ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio_AncFiltered.vcf

# Many filtering steps
# MQ or MQM = Mapping quality
# QUAL = Metric that is specific for the variant caller that denotes confidence
# DP = Read depth
# DP[x] or SAF,SAR,SRF,SRR = Array with read depth for Fwd, Rev, and from which strand
# Filters by quality, mapping quality, read depth, number of reads supporting variant, ballence between forward and reverse reads
(>2 echo ***BCFtools - Filter***)
bcftools filter -O v -o ${SAMPLE}_samtools_qualfiltered.vcf \
        -i 'MQ>30 & QUAL>75 & DP>10 & (DP4[2]+DP4[3])>4 & (DP4[2]+DP4[3])/DP>0.3 & (DP4[0]+DP4[2])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01 & (DP4[1]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01' \
        ${SAMPLE}_samtools_AB_AncFiltered.vcf

#Chris recommended a 40x coverage cutoff, but some of my samples are low coverage so I switched it to 10 for now.
bcftools filter -O v -o ${SAMPLE}_freebayes_qualfiltered.vcf \
        -i 'MQM>30 & MQMR>30 & QUAL>20 & INFO/DP>40 & (SAF+SAR)>4 & (SRF+SAF)/(INFO/DP)>0.01 & (SRR+SAR)/(INFO/DP)>0.01' \
        ${SAMPLE}_freebayes_BCBio_AncFiltered.vcf

# intersect samtools by freebayes to create final vcf to annotate
bedtools intersect -header \
        -a ${SAMPLE}_samtools_qualfiltered.vcf \
        -b ${SAMPLE}_freebayes_qualfiltered.vcf \
        > ${SAMPLE}_final.vcf

# Uses custom annotation script to put ORFs, tRNA, and ect. on the vcfs
(>2 echo ***Annotate***)
python3 ${SCRIPTS}/annotation_final.py \
        -f ${WORKDIR}/${SAMPLE}/${SAMPLE}_final.vcf \
        -s ${ANNOTATE}/orf_coding_all_R64-1-1_20110203.fasta \
        -n ${ANNOTATE}/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered \
        -g ${ANNOTATE}/S288C_reference_sequence_R64-1-1_20110203.fsa

# Remove intermediates
rm ${SAMPLE}_comb_R1R2.bam.intervals
rm ${SAMPLE}_comb_R1R2.MD.bam
rm ${SAMPLE}_comb_R1R2.RG.MD.bam
rm ${SAMPLE}_comb_R1R2.RG.MD.realign.bam
rm ${SAMPLE}_comb_R1R2.RG.MD.realign.bai
rm ${SAMPLE}_R1R2_sort.bam
rm ${SAMPLE}_R1R2_sort.bam.bai
rm ${SAMPLE}_comb_R1R2.RG.MD.sort.bam
rm ${SAMPLE}_comb_R1R2.RG.MD.sort.bam.bai

# remove intermediate samtools and freebayes files
rm ${SAMPLE}_samtools_AB_AncFiltered.vcf
rm ${SAMPLE}_freebayes_BCBio_AncFiltered.vcf

# remove some random files were somehow produced by the script (Not sure how it was produced)
# there is some file called 2 that was created. 
rm 2


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
module load zlib/1.2.13_compat
module load bwa/0.7.15
module load htslib/1.18
module load samtools/1.14
module load picard/3.1.1
module load GATK/3.7
module load python/3.12.1 numpy biopython lofreq/2.1.5-18
module load perl/5.26.3
module load VCFtools/0.1.16-20
module load bcftools/1.19
module load bedtools/2.25.0
module load freebayes/1.3.6

SAMPLE=$1 # Passed sample prefix (ex: Sample-01)
#ANC=$2
DIR=/net/dunham/vol2/Zilong/updating_pipeline_2024
WORKDIR=${DIR}/WorkDirectory # Where files will be created
SEQDIR=${DIR}/fastq # Location of Fastqs
SEQID=leah_FTevo # Project name and date for bam header
REF=/net/dunham/vol2/Zilong/updating_pipeline_2024/genomes/sacCer3.fasta # Reference genome
ANNOTATE=/net/dunham/vol2/Cris_L/ReferenceGenome/S288C_reference_genome_R64-1-1_20110203 # Location of custom annotation scripts
SCRIPTS=/net/dunham/vol2/Cris_L/Aaron_Reanalyze/Scripts # Location of custom scripts
#ANCBAM=${WORKDIR}/${ANC}/${ANC}_comb_R1R2.RG.MD.realign.sort.bam
#VCFDIR=${WORKDIR}/${ANC}/

cd ${WORKDIR}/${SAMPLE}/

# Many filtering steps
# MQ or MQM = Mapping quality
# QUAL = Metric that is specific for the variant caller that denotes confidence
# DP = Read depth
# DP[x] or SAF,SAR,SRF,SRR = Array with read depth for Fwd, Rev, and from which strand
# Filters by quality, mapping quality, read depth, number of reads supporting variant, ballence between forward and reverse reads
(>2 echo ***BCFtools - Filter***)
bcftools filter -O v -o ${SAMPLE}_samtools_filtered.vcf \
        -i 'MQ>30 & QUAL>75 & DP>10 & (DP4[2]+DP4[3])>4 & (DP4[2]+DP4[3])/DP>0.3 & (DP4[0]+DP4[2])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01 & (DP4[1]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01' \
        ${SAMPLE}_samtools_AB_AncFiltered.vcf

#module load numpy/1.7.0
#module load biopython/latest
#python ${DIR}/scripts/yeast_annotation_anna_edits.py 
#       -f ${SNPDIR}/${SAMPLE}_samtools_filtered.vcf 
#       -s ${DIR}/genomes/orf_coding_all_R64-1-1_20110203.fasta 
#       -n ${DIR}/genomes/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered 
#       -g ${DIR}/genomes/S288C_reference_sequence_R64-1-1_20110203.fsa


python ${SCRIPTS}/yeast_annotation_chris_edits_20170925.py \
        -f ${WORKDIR}/${SAMPLE}/${SAMPLE}_samtools_filtered.vcf \
        -s ${ANNOTATE}/orf_coding_all_R64-1-1_20110203.fasta \
        -n ${ANNOTATE}/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered \
        -g ${ANNOTATE}/S288C_reference_sequence_R64-1-1_20110203.fsa


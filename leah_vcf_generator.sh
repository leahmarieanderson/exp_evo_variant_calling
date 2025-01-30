#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Leah/fixing_pipeline_april2024
#$ -o /net/dunham/vol2/Leah/fixing_pipeline_april2024/outputs/
#$ -e /net/dunham/vol2/Leah/fixing_pipeline_april2024/errors/
#$ -l mfree=8G
#$ -l h_rt=36:0:0
#$ -N leah_vcf4_test_250129

## leah testing variant calling

module load modules modules-init modules-gs
module load zlib/1.3.1
module load bwa/0.7.17
module load htslib/1.19
module load samtools/1.19
module load picard/3.1.1
module load GATK/3.7
module load python/3.12.1 numpy biopython lofreq/2.1.5-18
module load perl/5.38.2
module load VCFtools/0.1.16-20
module load bcftools/1.20
module load bedtools/2.31.1
module load freebayes/1.3.6

FOLDER=fastq
SAMPLE=$1 # Passed sample prefix (ex: Sample-01)
ANC=$2
DIR=/net/dunham/vol2/Leah/fixing_pipeline_april2024
WORKDIR=${DIR}/WorkDirectory # Where files will be created
SEQDIR=${DIR}/${FOLDER} # Location of Fastqs
SEQID=leah_test # Project name and date for bam header
REF=${DIR}/genomes/sacCer3.fasta # Reference genome
ANNOTATE=${DIR}/genomes # Location of custom annotation scripts
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
ANCBAM=${WORKDIR}/${ANC}/bams/${ANC}_comb_R1R2.RG.MD.realign.sort.bam
VCFDIR=${WORKDIR}/${ANC}/

(>&2 echo ***LoFreq - Somatic***)
lofreq somatic -n ${ANCBAM} -t ${WORKDIR}/${SAMPLE}/bams/${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam -f ${REF} \
-o ${SAMPLE}_lofreq_

# Unzips lofreq vcfs
bgzip -d ${SAMPLE}_lofreq_somatic_final.snvs.vcf.gz
bgzip -d ${SAMPLE}_lofreq_tumor_relaxed.vcf.gz
bgzip -d ${SAMPLE}_lofreq_normal_relaxed.vcf.gz

# Filters samtools by ancestor
(>&2 echo ***Bedtools - Intersect***)
bedtools intersect -v -header \
        -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_samtools_AB.vcf \
        -b ${VCFDIR}/${ANC}_samtools_AB.vcf \
        > ${WORKDIR}/${SAMPLE}/${SAMPLE}_samtools_AB_AncFiltered.vcf

# Filters freebayes by ancestor 
bedtools intersect -v -header \
        -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio.vcf \
        -b ${VCFDIR}/${ANC}_freebayes_BCBio.vcf \
        > ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio_AncFiltered.vcf

# Filters lofreq by ancestor
bedtools intersect -v -header \
        -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_tumor_relaxed.vcf \
        -b ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_normal_relaxed.vcf \
        > ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_AncFiltered.vcf


# Annotate the AncFiltered
(>&2 echo ***Annotate***)
python3 ${SCRIPTS}/annotation_final.py \
        -f ${WORKDIR}/${SAMPLE}/${SAMPLE}_samtools_AB_AncFiltered.vcf \
        -s ${ANNOTATE}/orf_coding_all_R64-1-1_20110203.fasta \
        -n ${ANNOTATE}/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered \
        -g ${ANNOTATE}/S288C_reference_sequence_R64-1-1_20110203.fsa

python3 ${SCRIPTS}/annotation_final.py \
        -f ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio_AncFiltered.vcf \
        -s ${ANNOTATE}/orf_coding_all_R64-1-1_20110203.fasta \
        -n ${ANNOTATE}/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered \
        -g ${ANNOTATE}/S288C_reference_sequence_R64-1-1_20110203.fsa

python3 ${SCRIPTS}/annotation_final.py \
        -f ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_AncFiltered.vcf \
        -s ${ANNOTATE}/orf_coding_all_R64-1-1_20110203.fasta \
        -n ${ANNOTATE}/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered \
        -g ${ANNOTATE}/S288C_reference_sequence_R64-1-1_20110203.fsa


# Many filtering steps
# MQ or MQM = Mapping quality
# QUAL = Metric that is specific for the variant caller that denotes confidence
# DP = Read depth
# DP[x] or SAF,SAR,SRF,SRR = Array with read depth for Fwd, Rev, and from which strand
# Filters by quality, mapping quality, read depth, number of reads supporting variant, ballence between forward and reverse reads

(>&2 echo ***Apply Stringent Filter Based on Variant Caller and Return Combined CSV***)

# After, we would like to create a csv with just the necessary information
python3 ${SCRIPTS}/stringent_filter.py ${SAMPLE}_samtools_AB_AncFiltered_annotated_vcf.txt ${SAMPLE}_freebayes_BCBio_AncFiltered_annotated_vcf.txt ${SAMPLE}_lofreq_AncFiltered_annotated_vcf.txt

# remove all the lofreq intermediate files
rm ${SAMPLE}_lofreq_normal_relaxed.log
rm ${SAMPLE}_lofreq_normal_relaxed.vcf
rm ${SAMPLE}_lofreq_normal_relaxed.vcf.gz.tbi
rm ${SAMPLE}_lofreq_normal_stringent.indels.vcf.gz
rm ${SAMPLE}_lofreq_normal_stringent.indels.vcf.gz.tbi
rm ${SAMPLE}_lofreq_normal_stringent.snvs.vcf.gz
rm ${SAMPLE}_lofreq_normal_stringent.snvs.vcf.gz.tbi
rm ${SAMPLE}_lofreq_somatic_final.indels.vcf.gz
rm ${SAMPLE}_lofreq_somatic_final.indels.vcf.gz.tbi
rm ${SAMPLE}_lofreq_somatic_final.snvs.vcf
rm ${SAMPLE}_lofreq_somatic_final.snvs.vcf.gz.tbi
rm ${SAMPLE}_lofreq_somatic_raw.indels.vcf.gz
rm ${SAMPLE}_lofreq_somatic_raw.indels.vcf.gz.tbi
rm ${SAMPLE}_lofreq_somatic_raw.snvs.vcf.gz
rm ${SAMPLE}_lofreq_somatic_raw.snvs.vcf.gz.tbi
rm ${SAMPLE}_lofreq_tumor_relaxed.log
rm ${SAMPLE}_lofreq_tumor_relaxed.vcf
rm ${SAMPLE}_lofreq_tumor_relaxed.vcf.gz.tbi
rm ${SAMPLE}_lofreq_tumor_stringent.indels.vcf.gz
rm ${SAMPLE}_lofreq_tumor_stringent.indels.vcf.gz.tbi
rm ${SAMPLE}_lofreq_tumor_stringent.snvs.vcf.gz
rm ${SAMPLE}_lofreq_tumor_stringent.snvs.vcf.gz.tbi
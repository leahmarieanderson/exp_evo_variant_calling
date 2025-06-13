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
ANC=$2
DIR=/net/dunham/vol2/Zilong/updating_pipeline_2024
WORKDIR=${DIR}/WorkDirectory # Where files will be created
SEQDIR=${DIR}/${FOLDER} # Location of Fastqs
SCRIPTS=${DIR}/exp_evo_variant_calling # Path of annotation_final.py directory
SEQID=leah_freeze_evolution # Project name and date for bam header
REF=${SCRIPTS}/genomes/sacCer3.fasta # Reference genome
ANNOTATE=${SCRIPTS}/genomes # Location of custom annotation scripts
ANCBAM=${WORKDIR}/${ANC}/${ANC}_R1R2_MD.sort.bam
VCFDIR=${WORKDIR}/${ANC}

# In the case where ANC argument is given incorrectly, we don't want to take an hour creating files
# only to find out later on our code stopped working because we misspelled something or we didn't
# provide the correct file path for ANC

# Checks first if ANC name is provided in arguments
if [ -n "$2" ]; then 
        # Check if a file exists in the 'fastq' directory and contains ${ANC} in its name
        if find "$DIR/$FOLDER" -type f -name "${ANC}_*" | grep -q .; then
            (>&2 echo A fastq file with sample name ${ANC} exists in the ${FOLDER} directory.)
            # Check to see if ancestor.bam does not exists
            if [ ! -e ${ANCBAM} ]; then
            # Throw a warning that we cannot find the Ancestor bam and exit the shell script.
                (>&2 echo ***${ANCBAM} cannot be found***)
                (>&2 echo ***Need to create Ancestor bam or modify path***)
                exit 1 
            fi
        else
            (>&2 echo No fastq file with sample name ${ANC} was found in the ${FOLDER} directory.)
            (>&2 echo Possible misspelling of ancestor sample name or misplaced ancestor fastq file)
            (>&2 echo 'Check your fastq folder or double check that you did not misspell')
            # exit out 
            exit 1
        fi
fi

# Sets up folder structure
mkdir -p ${WORKDIR}/${SAMPLE}
cd ${WORKDIR}/${SAMPLE}

# Perform FastQC checks on our samples 
(>&2 echo ***FASTQC on SAMPLE ***)
fastqc ${SEQDIR}/${SAMPLE}_*R1*.fastq.gz -o ${WORKDIR}/${SAMPLE}/
fastqc ${SEQDIR}/${SAMPLE}_*R2*.fastq.gz -o ${WORKDIR}/${SAMPLE}/

# remove the zip file since the html file has everything we need to know
rm ${SAMPLE}_*R1*fastqc.zip
rm ${SAMPLE}_*R2*fastqc.zip

# Align reads with bwa
(>&2 echo ***BWA - mem -R***)
bwa mem -R "@RG\tID:${SEQID}\tSM:${SAMPLE}\tLB:1" ${REF} ${SEQDIR}/${SAMPLE}_*R1*.fastq.gz ${SEQDIR}/${SAMPLE}_*R2*.fastq.gz > ${SAMPLE}_R1R2.sam

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

# for now, we don't know the known sites for variants for our yeast, so we will skip this step
# (>&2 echo ***GATK4 - BaseRecalibrator***)

# 1. build the model
#gatk BaseRecalibrator \
#        -I ${SAMPLE}_R1R2_MD.sort.bam \
#        -R ${REF} --known-sites ${known_sites} \
#        -O ${SAMPLE}_recal_data.table


# 2. Apply the model to adjust the base quality scores
#gatk ApplyBQSR \
#        -I ${SAMPLE}_R1R2_MD.sort.bam \
#        -R ${REF} 
#        --bqsr-recal-file ${SAMPLE}_recal_data.table 
#        -O ${SAMPLE}_R1R2_MD.sort.bqrs.bam

(>&2 echo ***GATK4 - Calling Variants***)
gatk HaplotypeCaller \
     -R ${REF} \
     -I ${SAMPLE}_R1R2_MD.sort.bam \
     -O ${SAMPLE}_gatk_haplo.vcf

# Freebayes with a lot of arguments for population calling
freebayes -f ${REF} \
        --pooled-continuous --report-genotype-likelihood-max --allele-balance-priors-off --min-alternate-fraction 0.1 \
        ${SAMPLE}_R1R2_MD.sort.bam > ${SAMPLE}_freebayes_BCBio.vcf

# Requires ANC from this line down
# check if ANC argument was given. If there is, then continue with Ancestor filtering 
if [ -n "$2" ]; then 
        # Go to Work Directory
        cd ${WORKDIR}/${SAMPLE}
        (>&2 echo ***LoFreq - Somatic***)
        lofreq somatic -n ${ANCBAM} -t ${WORKDIR}/${SAMPLE}/${SAMPLE}_R1R2_MD.sort.bam -f ${REF} \
        -o ${SAMPLE}_lofreq_

        # Unzips lofreq vcfs
        bgzip -d ${SAMPLE}_lofreq_somatic_final.snvs.vcf.gz
        bgzip -d ${SAMPLE}_lofreq_tumor_relaxed.vcf.gz
        bgzip -d ${SAMPLE}_lofreq_normal_relaxed.vcf.gz

        # Filters gatk_haplo by ancestor
        (>&2 echo ***Bedtools - Intersect***)
        bedtools intersect -v -header \
                -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_gatk_haplo.vcf \
                -b ${VCFDIR}/${ANC}_gatk_haplo_quality_filter.vcf \
                > ${WORKDIR}/${SAMPLE}/${SAMPLE}_gatk_haplo_AncFiltered.vcf

        # Filters freebayes by ancestor 
        bedtools intersect -v -header \
                -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio.vcf \
                -b ${VCFDIR}/${ANC}_freebayes_BCBio_quality_filter.vcf \
                > ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio_AncFiltered.vcf

        # Filters lofreq by ancestor
        bedtools intersect -v -header \
                -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_tumor_relaxed.vcf \
                -b ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_normal_relaxed.vcf \
                > ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_AncFiltered.vcf

        # Annotate the AncFiltered
        (>&2 echo ***Annotate***)
        python3 ${SCRIPTS}/annotation_final.py \
                -f ${WORKDIR}/${SAMPLE}/${SAMPLE}_gatk_haplo_AncFiltered.vcf \
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
        # SOR - Symmetric Odds Ratio for strain bias, lower values are better
        # or SAF,SAR,SRF,SRR = Array with read depth for Fwd, Rev, and from which strand
        # Filters by quality, mapping quality, read depth, number of reads supporting variant, ballence between forward and reverse reads

        (>&2 echo ***Apply Stringent Filter Based on Variant Caller and Return Combined CSV***)

        # After, we would like to create a csv with just the necessary information
        python3 ${SCRIPTS}/stringent_filter.py \
        ${SAMPLE}_gatk_haplo_AncFiltered_annotated_vcf.txt \
        ${SAMPLE}_freebayes_BCBio_AncFiltered_annotated_vcf.txt \
        ${SAMPLE}_lofreq_AncFiltered_annotated_vcf.txt

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

else 
# use a quality and read depth filter on the ancestor vcfs

(>&2 echo ***BCFtools - Filter***)
bcftools filter -O v -o ${SAMPLE}_gatk_haplo_quality_filter.vcf \
        -i 'MQ>30 & QUAL>75 & (INFO/DP)>40 & (INFO/SOR)<3' \
        ${SAMPLE}_gatk_haplo.vcf

bcftools filter -O v -o ${SAMPLE}_freebayes_BCBio_quality_filter.vcf \
        -i 'MQM>30 & QUAL>20 & INFO/DP>10 & (SAF+SAR)>4 & (SRF+SAF)/(INFO/DP)>0.01 & (SRR+SAR)/(INFO/DP)>0.01' \
        ${SAMPLE}_freebayes_BCBio.vcf
fi

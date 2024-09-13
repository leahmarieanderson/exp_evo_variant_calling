#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Leah/230331_FTevo
#$ -e /net/dunham/vol2/Leah/230331_FTevo/errors/
#$ -o /net/dunham/vol2/Leah/230331_FTevo/outputs/
#$ -N LMA_H15
#$ -l mfree=4G

## CNV pipeline for figuring out CN from bam files: uses both wig file from igvtools and mpileup

module load modules modules-init modules-gs
module load python/2.7.13
module load java/1.8.0
module load GATK/3.7
IGVTOOLS=/net/dunham/vol2/Caiti/hybrid_seq/IGVTools/igvtools.jar

## sh runCNScript.sh SAMPLE 1000 ANCESTOR 1
#cer

SAMPLE=$1 #sample prefix (ex: Sample-01)
SIZE=$2
WORKDIR=/net/dunham/vol2/Leah/230331_FTevo/WorkDirectory
BAMDIR=${WORKDIR}/${SAMPLE}
CNDIR=/net/dunham/vol2/Leah/230331_FTevo/WorkDirectory  #CHANGE back to ${WORKDIR}/${SAMPLE}/CNV_new_${SIZE}bp
SCRIPTS=/net/dunham/vol2/Bryce/sequencing_data/Rowley_Project/Scripts
REF=/net/dunham/vol2/Caiti/reference_seq/sacCer3.fasta
ANC=$3
PLOIDY=$4  #change back to $4 if have ancestor
CERORF=/net/dunham/vol2/Caiti/reference_seq/cer_homolog_coordinates_filt.txt

cd ${WORKDIR}/${SAMPLE}
mkdir -p CNV_new_${SIZE}bp

## Get depth of coverage info (old version ran this on earlier bam). Ignores mito
java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar -T DepthOfCoverage \
	-R ${REF} -I ${BAMDIR}/${SAMPLE}_comb_R1R2.RG.MD.sort.bam \
	-o ${CNDIR}/${SAMPLE}_comb_R1R2.RG.MD.sort.bam.DOC \
	-XL chrM -omitBaseOutput -omitLocusTable -omitIntervals -rf BadCigar

## Make wig file (can also run igvtools directly, but java allows mem management)
## Can change window size
## NOTE: RUN ON ANCESTRAL FIRST
java -Xmx2g -Djava.awt.headless=true -jar $IGVTOOLS count -w ${SIZE} --minMapQuality 30 \
	${BAMDIR}/${SAMPLE}_comb_R1R2.RG.MD.sort.bam \
	${CNDIR}/${SAMPLE}_${SIZE}bp.wig ${REF}
java -Xmx2g -Djava.awt.headless=true -jar $IGVTOOLS count -w ${SIZE} --minMapQuality 0 \
        ${BAMDIR}/${SAMPLE}_comb_R1R2.RG.MD.sort.bam \
        ${CNDIR}/${SAMPLE}_${SIZE}bp.All.wig ${REF}
python ${SCRIPTS}/wigNormalizedToAverageReadDepth_MapQ_ForPlot.py \
        ${CNDIR}/${SAMPLE}_comb_R1R2.RG.MD.sort.bam.DOC.sample_summary \
        ${CNDIR}/${SAMPLE}_${SIZE}bp.wig \
	${CNDIR}/${SAMPLE}_${SIZE}bp.All.wig \
	${PLOIDY} \
        ${CNDIR}/${SAMPLE}_${SIZE}bp_norm.wig

module load gcc/8.1.0
module load R/latest

Rscript ${SCRIPTS}/PlotCopyNumber_OneSample_20200106.R ${WORKDIR} ${SAMPLE}
Rscript ${SCRIPTS}/PlotCopyNumber_TwoSample_20200106.R ${WORKDIR} ${SAMPLE} ${ANC}

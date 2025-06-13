#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Leah/labmeeting_250613
#$ -o /net/dunham/vol2/Leah/labmeeting_250613/outputs/
#$ -e /net/dunham/vol2/Leah/labmeeting_250613/errors/
#$ -N labMeeting1
#$ -l mfree=4G

## CNV pipeline for figuring out CN from bam files: uses both wig file from igvtools and mpileup

module load modules modules-init modules-gs
module load python/3.12.1
module load java/1.8.0
module load GATK/3.7
IGVTOOLS=/net/dunham/vol2/Caiti/hybrid_seq/IGVTools/igvtools.jar

## sh runCNScript.sh SAMPLE 1000 1 ANCESTOR if you want the two sample plot
## else just do 
## sh runCNScript.sh SAMPLE 1000 1
#cer

SAMPLE=$1 #sample prefix (ex: Sample-01)
SIZE=$2
DIR=/net/dunham/vol2/Leah/labmeeting_250613
WORKDIR=${DIR}/WorkDirectory
BAMDIR=${WORKDIR}/${SAMPLE}
CNDIR=${WORKDIR}/${SAMPLE}/CNV_new_${SIZE}bp  #CHANGE back to ${WORKDIR}/${SAMPLE}/CNV_new_${SIZE}bp
SCRIPTS=${DIR}/exp_evo_variant_calling
REF=/net/dunham/vol2/Caiti/reference_seq/sacCer3.fasta
PLOIDY=$3  
ANC=$4
# CERORF=/net/dunham/vol2/Caiti/reference_seq/cer_homolog_coordinates_filt.txt # we aren't using this at the moment

cd ${WORKDIR}/${SAMPLE}
mkdir -p CNV_new_${SIZE}bp

## Get depth of coverage info (old version ran this on earlier bam). Ignores mito
java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar -T DepthOfCoverage \
	-R ${REF} -I ${BAMDIR}/${SAMPLE}_R1R2_MD.sort.bam \
	-o ${CNDIR}/${SAMPLE}_R1R2_MD.sort.bam.DOC \
	-XL chrM -omitBaseOutput -omitLocusTable -omitIntervals -rf BadCigar

## Make wig file (can also run igvtools directly, but java allows mem management)
## Can change window size
# Make a wig file with data that satisfies the minimum mapping quality 
java -Xmx2g -Djava.awt.headless=true -jar $IGVTOOLS count -w ${SIZE} --minMapQuality 30 \
	${BAMDIR}/${SAMPLE}_R1R2_MD.sort.bam \
	${CNDIR}/${SAMPLE}_${SIZE}bp.wig ${REF}

# Make the wig file that contains everything
java -Xmx2g -Djava.awt.headless=true -jar $IGVTOOLS count -w ${SIZE} --minMapQuality 0 \
    ${BAMDIR}/${SAMPLE}_R1R2_MD.sort.bam \
    ${CNDIR}/${SAMPLE}_${SIZE}bp.All.wig ${REF}
python ${SCRIPTS}/wigNormalizedToAverageReadDepth_MapQ_ForPlot.py \
    ${CNDIR}/${SAMPLE}_R1R2_MD.sort.bam.DOC.sample_summary \
    ${CNDIR}/${SAMPLE}_${SIZE}bp.wig \
	${CNDIR}/${SAMPLE}_${SIZE}bp.All.wig \
	${PLOIDY} \
    ${CNDIR}/${SAMPLE}_${SIZE}bp_norm.wig

module load R/latest

# if there isn't an argument for ANCESTOR
if [ -z "$4" ]; then 
	# There isn't an argument for the ancestor, just plot our sample
	cd ${WORKDIR}/${SAMPLE}
	Rscript ${SCRIPTS}/PlotCopyNumber_OneSample_20200106.R ${WORKDIR} ${SAMPLE} ${SIZE}
else
	# There is an argument for the ancestor
	# check directories to see if it has the .wig file needed for PlotCopyNumber_TwoSample_20200106.R
	file_path=${WORKDIR}/${ANC}/CNV_new_${SIZE}bp/${ANC}_${SIZE}bp_norm.wig
	if [ ! -f $file_path ]; then # if ancestor bp_norm.wig file does not exists then create it
		# create the directories and files necessary to run the two_sample R script.
		# basically the same thing we did for our sample but now with the Ancestor 
		cd ${WORKDIR}/${ANC}
		mkdir -p CNV_new_${SIZE}bp

		# Creating Ancestor Wig files

		java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar -T DepthOfCoverage \
			-R ${REF} -I ${WORKDIR}/${ANC}/${ANC}_R1R2_MD.sort.bam \
			-o ${WORKDIR}/${ANC}/CNV_new_${SIZE}bp/${ANC}_R1R2_MD.sort.bam.DOC \
			-XL chrM -omitBaseOutput -omitLocusTable -omitIntervals -rf BadCigar

		java -Xmx2g -Djava.awt.headless=true -jar $IGVTOOLS count -w ${SIZE} --minMapQuality 30 \
			${WORKDIR}/${ANC}/${ANC}_R1R2_MD.sort.bam \
			${WORKDIR}/${ANC}/CNV_new_${SIZE}bp/${ANC}_${SIZE}bp.wig ${REF}

		java -Xmx2g -Djava.awt.headless=true -jar $IGVTOOLS count -w ${SIZE} --minMapQuality 0 \
			${WORKDIR}/${ANC}/${ANC}_R1R2_MD.sort.bam \
			${WORKDIR}/${ANC}/CNV_new_${SIZE}bp/${ANC}_${SIZE}bp.All.wig ${REF}

		python ${SCRIPTS}/wigNormalizedToAverageReadDepth_MapQ_ForPlot.py \
    		${WORKDIR}/${ANC}/CNV_new_${SIZE}bp/${ANC}_R1R2_MD.sort.bam.DOC.sample_summary \
    		${WORKDIR}/${ANC}/CNV_new_${SIZE}bp/${ANC}_${SIZE}bp.wig \
			${WORKDIR}/${ANC}/CNV_new_${SIZE}bp/${ANC}_${SIZE}bp.All.wig \
			${PLOIDY} \
    		${WORKDIR}/${ANC}/CNV_new_${SIZE}bp/${ANC}_${SIZE}bp_norm.wig

		# After creating all required files, now plot both sample then sample with ancestor. 
		cd ${WORKDIR}/${SAMPLE}
		Rscript ${SCRIPTS}/PlotCopyNumber_OneSample_20200106.R ${WORKDIR} ${SAMPLE} ${SIZE}
		Rscript ${SCRIPTS}/PlotCopyNumber_TwoSample_20200106.R ${WORKDIR} ${SAMPLE} ${SIZE} ${ANC} 

	else 
		cd ${WORKDIR}/${SAMPLE}
		Rscript ${SCRIPTS}/PlotCopyNumber_OneSample_20200106.R ${WORKDIR} ${SAMPLE} ${SIZE}
		cd ${WORKDIR}/${SAMPLE}
		Rscript ${SCRIPTS}/PlotCopyNumber_TwoSample_20200106.R ${WORKDIR} ${SAMPLE} ${SIZE} ${ANC} 
	fi
fi


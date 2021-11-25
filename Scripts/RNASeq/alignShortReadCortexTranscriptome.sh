## process short read data and align to CORTEX TRANSCRIPTOME


sampleName=$1
outPrefix=CortexMergedSMRT

cd ${RNASeqDIR}

TRIMFOLDER=Trimmed

star_f1=$(ls ${TRIMFOLDER}/${sampleName}_*[rR]1*f*z)
star_f2=$(ls ${TRIMFOLDER}/${sampleName}_*[rR]2*f*z)

mkdir -p ${ALIGNEDRNA}

module purge
module load STAR

if [ ! -f ${ALIGNEDRNA}/${outPrefix}_${sampleName}Log.final.out ]
then	
	echo "aligning" ${star_f1}

	## align with STAR using GENCODE v38 star index
	STAR --genomeDir ${STARINDEXDIR}/Brain --runThreadN 18 --readFilesIn ${star_f1},${star_f2} \
		--readFilesCommand zcat \
		--outFileNamePrefix ${ALIGNEDRNA}/${outPrefix}_${sampleName} \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes Standard
fi


	
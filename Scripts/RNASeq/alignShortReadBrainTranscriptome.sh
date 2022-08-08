
sampleName=$1
RNASEQDIR=$2
ALIGNEDDIR=$3
GENECOUNTDIR=$4
TRANSCRIPTOME=$5

TRIMDIR=trimmed

##find trimmed files

cd ${RNASEQDIR}
star_f1=$(ls ${TRIMDIR}/${sampleName}*[rR]1*f*z)
star_f2=$(ls ${TRIMDIR}/${sampleName}*[rR]2*f*z)

mkdir -p ${ALIGNEDDIR}/${TRANSCRIPTOME}


if [ ! -f ${ALIGNEDDIR}/${TRANSCRIPTOME}/${TRANSCRIPTOME}_${sampleName}Log.final.out ]
then	
	echo "aligning" ${star_f1}

	## align with STAR to transcriptome provided on command line
	STAR --genomeDir ${STARINDEXDIR}/${TRANSCRIPTOME} --runThreadN 18 --readFilesIn ${star_f1},${star_f2} \
		--readFilesCommand zcat \
		--outFileNamePrefix ${ALIGNEDDIR}/${TRANSCRIPTOME}/${TRANSCRIPTOME}_${sampleName} \
		--outSAMtype BAM Unsorted \
		--outSAMunmapped Within \
		--outSAMattributes  NH HI AS NM MD \
		--outFilterType BySJout \
		--outFilterMultimapNmax 20 \
		--outFilterMismatchNmax 999  \
		--outFilterMismatchNoverLmax 0.04  \
		--alignIntronMin 20  \
		--alignIntronMax 1000000  \
		--alignMatesGapMax 1000000  \
		--alignSJoverhangMin 8  \
		--alignSJDBoverhangMin 1  \
		--sjdbScore 1  \
		--genomeLoad NoSharedMemory  \
		--quantMode TranscriptomeSAM  \
		--outSAMheaderHD \@HD VN:1.4 SO:unsorted \
		--limitOutSJcollapsed 2000000 

fi

mkdir -p ${GENECOUNTDIR}/RSEM/${TRANSCRIPTOME}

if [ ! -f ${GENECOUNTDIR}/RSEM/${TRANSCRIPTOME}/${sampleName}.isoforms.results ]
then
	rsem-calculate-expression --num-threads 10 --alignments ${ALIGNEDDIR}/${TRANSCRIPTOME}/${TRANSCRIPTOME}_${sampleName}Aligned.toTranscriptome.out.bam ${RSEMREFDIR}/${TRANSCRIPTOME}/${TRANSCRIPTOME} ${GENECOUNTDIR}/RSEM/${TRANSCRIPTOME}/${sampleName}
	rm ${ALIGNEDDIR}/${TRANSCRIPTOME}/${TRANSCRIPTOME}_${sampleName}Aligned.out.bam
fi
	
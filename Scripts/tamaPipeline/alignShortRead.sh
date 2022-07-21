## SQANTI requires short read data aligned to the isoseq transcriptome
## this script generates gene/isoform counts matrix using RSEM which internally generates the STAR index for personalized transcriptome and aligns RNA seq data to it 
## had to split into separate STAR and RSEM steps as STAR default --limitOutSJcollapsed parameter was too low.
## recommended to run each sample in turn rather than as a single sample
## find trimmed files


star_f1=$1
star_f2=$2
rnaID=$(basename ${star_f1})
rnaID=${rnaID%%_*}

echo "Aligning short read data for sample " ${rnaID}

if [ ! -f ${GENECOUNTPATH}/RSEM/TAMAmerge/${rnaID}.isoforms.results ]
then

	STAR --genomeDir ${RSEMREFDIR}/TAMAmerge/ --runThreadN 18 --readFilesIn ${star_f1},${star_f2} \
		--readFilesCommand zcat \
		--outFileNamePrefix ${ALIGNEDDIR}/TAMAmerge/${rnaID} \
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

	rsem-calculate-expression --num-threads 10 --alignments ${ALIGNEDDIR}/TAMAmerge/${rnaID}Aligned.toTranscriptome.out.bam ${RSEMREFDIR}/TAMAmerge/TAMAmerge ${GENECOUNTPATH}/RSEM/TAMAmerge/${rnaID}
	rm ${ALIGNEDDIR}/TAMAmerge/${rnaID}Aligned.out.bam
fi



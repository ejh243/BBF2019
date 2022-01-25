## SQANTI requires short read data aligned to the isoseq transcriptome
## this script generates gene/isoform counts matrix using RSEM which internally generates the STAR index for personalized transcriptome and aligns RNA seq data to it 
## had to split into separate STAR and RSEM steps as STAR default --limitOutSJcollapsed parameter was too low.

p=$1
basename=${p%.subreads.bam}
echo "Aligning short read data to " ${basename} "transcriptome"

## find trimmed files store incomma separated list
star_f1=$(ls -m ${RNASeqDIR}Trimmed/*[rR]1*f*z | tr -d '\n')
star_f2=$(ls -m ${RNASeqDIR}Trimmed/*[rR]2*f*z | tr -d '\n')

GFF=${ALIGNEDDIR}/Collapsed/${basename}/out.collapsed.filtered.gff

## use RSEM to generate counts matrix

if [ ! -f ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/${basename}.isoforms.results ]
then
    mkdir -p ${RSEMREFDIR}/${basename}/
    rsem-prepare-reference --gtf ${GFF} --star ${REFGENOME} ${RSEMREFDIR}/${basename}/${basename}
    
    #rsem-calculate-expression --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/${basename}/${basename} ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/${basename}
    
    ## align with STAR 
    STAR --genomeDir ${RSEMREFDIR}/${basename}/ --runThreadN 18 --readFilesIn ${star_f1},${star_f2} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${RSEMREFDIR}/${basename}/${basename} \
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
		
	rm ${RSEMREFDIR}/${basename}/${basename}Aligned.out.bam	
    rsem-calculate-expression --num-threads 10 --alignments ${RSEMREFDIR}/${basename}/${basename}Aligned.toTranscriptome.out.bam ${RSEMREFDIR}/${basename}/${basename} ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/${basename}
	
	
fi



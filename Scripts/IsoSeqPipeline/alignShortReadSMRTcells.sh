## SQANTI requires short read data aligned to the isoseq transcriptome
## this script generates gene/isoform counts matrix using RSEM which internally generates the STAR index for personalized transcriptome and aligns RNA seq data to it 

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
	
	rsem-calculate-expression --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/${basename}/${basename} ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/${basename}
fi


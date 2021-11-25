#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/alignPersonalTranscriptome
#SBATCH --error=LogFiles/alignPersonalTranscriptome
#SBATCH --job-name=alignPersonalTranscriptome

## to run SQANTI requires short read data from that sample aligned to ISOSeq transcriptome

## this script generates gene/isoform counts matrix using RSEM which internally generates the STAR index for personalized transcriptome and aligns RNA seq data to it 

## take SMRT cell level isoforms


source ./Config/config.txt
module load STAR
module load RSEM

mkdir -p ${RSEMREFDIR}
mkdir -p ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/

## load matched sample ids from sample sheet
{
read ## to ignore first line
while IFS=',' read -ra ADDR
do
  isoseqID=${ADDR[2]%.subreads.bam}
  rnaID=${ADDR[8]}
  
  echo "Iso-Seq Sample" $isoseqID "is matched with RNA-Seq Sample" $rnaID 

	##find trimmed files
	star_f1=$(ls ${RNASeqDIR}Trimmed/${rnaID}*[rR]1*f*z)
	star_f2=$(ls ${RNASeqDIR}Trimmed/${rnaID}*[rR]2*f*z)

	GFF=${ALIGNEDDIR}/Collapsed/${isoseqID}/out.collapsed.filtered.gff

	## use RSEM to generate counts matrix

	if [ ! -f ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/${isoseqID}.isoforms.results ]
	then
		mkdir -p ${RSEMREFDIR}/${isoseqID}/${isoseqID}
		rsem-prepare-reference --gtf ${GFF} --star ${REFGENOME} ${RSEMREFDIR}/${isoseqID}/${isoseqID}
	
		rsem-calculate-expression --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/${isoseqID}/${isoseqID} ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/${isoseqID}
	fi

done } < "$SAMPLESHEET"
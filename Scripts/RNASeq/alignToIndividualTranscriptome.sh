#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/alignIndividualTranscriptome
#SBATCH --error=LogFiles/alignIndividualTranscriptome
#SBATCH --job-name=alignIndividualTranscriptome

## to run SQANTI requires short read data from that sample aligned to ISOSeq transcriptome

## this script generates gene/isoform counts matrix using RSEM which internally generates the STAR index for personalized transcriptome and aligns RNA seq data to it 

## takes individual cell level isoforms (i.e. where multiple SMRT cells merged together before clustering)


source ./Config/config.txt
module load STAR
module load RSEM

GFFFILES=$(ls ${ALIGNEDDIR}/Collapsed/[0-9]*_*/out.collapsed.filtered.gff)

for GFF in ${GFFFILES[@]}
do
	FOLDER=$(basename ${GFF%/out.collapsed.filtered.gff})
	rnaID=${FOLDER%%_*}
	
	echo ${rnaID}
	
	##find trimmed files
	star_f1=$(ls ${RNASeqDIR}Trimmed/${rnaID}*[rR]1*f*z)
	star_f2=$(ls ${RNASeqDIR}Trimmed/${rnaID}*[rR]2*f*z)
	
	## use RSEM to generate counts matrix
	mkdir -p ${RSEMREFDIR}/${FOLDER}/${FOLDER}

	rsem-prepare-reference --gtf ${GFF} --star ${REFGENOME} ${RSEMREFDIR}/${FOLDER}/${FOLDER}

	mkdir -p ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/
	rsem-calculate-expression -p 8 --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/${FOLDER}/${FOLDER} ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/${FOLDER}

done

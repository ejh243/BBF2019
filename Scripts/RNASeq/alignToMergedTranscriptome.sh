#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/alignMergedTranscriptome
#SBATCH --error=LogFiles/alignMergedTranscriptome
#SBATCH --job-name=alignMergedTranscriptome

## this script generates gene/isoform counts matrix using RSEM which internally generates the STAR index for personalized transcriptome and aligns RNA seq data to it 

## take merged, polished isoforms
## align all short read data


source ./Config/config.txt
module load STAR
module load RSEM

GFF=${ALIGNEDDIR}/Collapsed/merged.collapsed.filtered.gff
	
##find trimmed files
star_f1=($(ls -m ${RNASeqDIR}Trimmed/*[rR]1*f*z))
star_f2=($(ls -m ${RNASeqDIR}Trimmed/*[rR]2*f*z))

star_f1=$(printf "%s" "${star_f1[@]}" && echo "")
star_f2=$(printf "%s" "${star_f2[@]}" && echo "")

## use RSEM to generate counts matrix
mkdir -p ${RSEMREFDIR}

#rsem-prepare-reference --gtf ${GFF} --star ${REFGENOME} ${RSEMREFDIR}/merged

mkdir -p ${GENECOUNTSDIR}/RSEM/MergedTranscriptome/
rsem-calculate-expression -p 8 --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/merged ${GENECOUNTSDIR}/RSEM/MergedTranscriptome/All


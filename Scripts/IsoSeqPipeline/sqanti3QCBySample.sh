#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the high memory queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/sqanti.o
#SBATCH --error=LogFiles/sqanti.e
#SBATCH --job-name=sqanti

## prior to merging to get master transcriptome use sqanti3 to qc sample-level transcriptomes

## require short read data aligned to individual transcriptomes

module load Miniconda2
source ./Config/config.txt
source activate SQANTI3.env

## input is collapsed high quality isoforms
## this are organized into different folders
## as this script is 
sampleDIR=($(ls -d ${ALIGNEDDIR}/Collapsed/*/))

for sample in ${sampleDIR[@]}
do
	sampleName=$(basename ${sample})
	
	## check if gene expression matrix for this transcriptome exists
	if [ -f ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/${sampleName}.isoforms.results ]
	then
		echo "Running SQANTI3 on" ${sampleName}
	
		mkdir -p ${sample}/SQANTI3
		python ${SOFTWAREPATH}/SQANTI3/sqanti3_qc.py ${sample}/out.collapsed.filtered.gff ${GENCODEGTF} ${REFGENOME} \
		-o $sampleName -d ${sample}/SQANTI3/ --cage_peak ${SQANTIDATA}/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed    \
			--polyA_motif_list ${SQANTIDATA}/polyA_motifs/mouse_and_human.polyA_motif.txt    \
			 -fl ${sample}/out.collapsed.filtered.abundance.txt  --isoAnnotLite --saturation --expression ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/${sampleName}.isoforms.results --cpus 4 --report both 
	else
		echo "Expression counts matrix does not exist for " ${sampleName} ": skipping"
	fi

done
	


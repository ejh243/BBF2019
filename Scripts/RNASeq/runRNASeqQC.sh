#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the high memory queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/rnaseqc.o
#SBATCH --error=LogFiles/rnaseqc.e
#SBATCH --job-name=rnaseqc

module load Miniconda2
source ./Config/config.txt
source activate rnaseqc

## create clustered gencode gtf for 
## using script from https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model
python ${SOFTWAREPATH}gtex-pipeline/collapse_annotation.py ${GENCODEGTF} ${GENCODEGTF/%annotation.gtf/genes.gtf}

cd ${SOFTWAREPATH}/rnaseqc

## run rnaseqc for each bam file
mkdir -p ${ALIGNEDRNA}/rnaseqc
for bamfile in `ls ${ALIGNEDRNA}/*.bam`; 
do
	sampleName=$(basename ${bamfile})
	sampleName=${sampleName%Aligned.sortedByCoord.out.bam}
	sampleName=${sampleName#GENCODE_v38_}
	 ${SOFTWAREPATH}/rnaseqc.v2.4.2.linux ${GENCODEGTF/%annotation.gtf/genes.gtf} ${bamfile} ${ALIGNEDRNA}/rnaseqc/ -s ${sampleName}
done

## aggregate results
mkdir -p ${ALIGNEDRNA}/rnaseqc/merged
python3 -m rnaseqc aggregate -o ${ALIGNEDRNA}/rnaseqc/merged ${ALIGNEDRNA}/rnaseqc/ LBB

## create figures
python3 -m rnaseqc report --tpm ${ALIGNEDRNA}/rnaseqc/merged/LBB.gene_tpm.gct.gz  --output-dir ${ALIGNEDRNA}/rnaseqc/merged ${ALIGNEDRNA}/rnaseqc/merged/LBB.metrics.txt.gz LBB


#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the high memory queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=LogFiles/CMCprocess-%A_%a.o
#SBATCH --error=LogFiles/CMCprocess-%A_%a.e
#SBATCH --job-name=CMCprocess-%A_%a

## this script converts bam file into fastq files
INDIR=$1

mkdir -p ${INDIR}/fastq_tmp/
mkdir -p ${INDIR}/fastq
BAMFILES=($(ls ${INDIR}*.bam))

bam_file=${BAMFILES[${SLURM_ARRAY_TASK_ID}]}

sampleName=$(basename ${bam_file%.accepted_hits.sort.coord.bam})

module load SAMtools

samtools collate $bam_file ${INDIR}/fastq_tmp/${sampleName}.collated
samtools fastq -F 2816 -c 6 -1 ${INDIR}/fastq_tmp/${sampleName}_1.fastq.gz -2 ${INDIR}/fastq_tmp/${sampleName}_2.fastq.gz ${INDIR}/fastq_tmp/${sampleName}.collated.bam

## merge with unmapped reads
cat ${INDIR}/fastq_tmp/${sampleName}_1.fastq.gz ${INDIR}/${sampleName}_accepted_hitsUnmapped.out.mate1.gz > ${INDIR}/fastq/${sampleName}_merged_1.fastq.gz
cat ${INDIR}/fastq_tmp/${sampleName}_2.fastq.gz ${INDIR}/${sampleName}_accepted_hitsUnmapped.out.mate2.gz > ${INDIR}/fastq/${sampleName}_merged_2.fastq.gz

rm ${INDIR}/fastq_tmp/${sampleName}*


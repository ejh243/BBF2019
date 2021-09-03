#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p sq # submit to the high memory queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/AlignFilter.o
#SBATCH --error=LogFiles/AlignFilter.e
#SBATCH --job-name=MergeAndCluster

## this script processes data where all samples have been combined see alternative script to process each sample individually

cd $SLURM_SUBMIT_DIR

module load Miniconda2
source ./Config/config.txt


echo "Changing Folder to: "
echo $DATADIR

module load minimap2

## align with minimap2
minimap2 -ax splice -t 30 -uf --secondary=no -C5 \
   ${REFGENOME} ${PROCESSEDDIR}/Cluster/clustered.hq.fasta.gz > ${ALIGNEDDIR}/clustered.hq_isoforms.fastq.sam

sort -k 3,3 -k 4,4n hq_isoforms.fastq.sam > hq_isoforms.fastq.sorted.sam

## collapse redundant isoforms
mkdir -p ${ALIGNEDDIR}/Collapsed/

source activate anaCogent

collapse_isoforms_by_sam.py --input ${PROCESSEDDIR}/Cluster/clustered.hq.fasta.gz --fq \
   -s ${ALIGNEDDIR}/clustered.hq_isoforms.fastq.sam --dun-merge-5-shorter -o ${ALIGNEDDIR}/Collapsed/clustered
   
## generate counts

get_abundance_post_collapse.py ${ALIGNEDDIR}/Collapsed/clustered.collapsed ${PROCESSEDDIR}/Cluster/clustered.cluster_report.csv
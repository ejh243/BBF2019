#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the high memory queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/AlignFilter.o
#SBATCH --error=LogFiles/AlignFilter.e
#SBATCH --job-name=AlignFilter

## this script processes data where all samples have been combined see alternative script to process each sample individually

cd $SLURM_SUBMIT_DIR


source ./Config/config.txt


echo "Changing Folder to: "
echo $DATADIR

module load minimap2

## align with minimap2
minimap2 -ax splice -t 30 -uf --secondary=no -C5 \
   ${REFGENOME} ${PROCESSEDDIR}/Polish/polished.hq.fastq.gz > ${ALIGNEDDIR}/polished.hq_isoforms.sam

sort -k 3,3 -k 4,4n ${ALIGNEDDIR}/polished.hq_isoforms.sam > ${ALIGNEDDIR}/polished.hq_isoforms.sorted.sam

## collapse redundant isoforms
mkdir -p ${ALIGNEDDIR}/Collapsed/


module load Miniconda2
source activate anaCogent

gunzip ${PROCESSEDDIR}/Polish/polished.hq.fastq.gz

collapse_isoforms_by_sam.py --input ${PROCESSEDDIR}/Polish/polished.hq.fastq --fq \
   -s ${ALIGNEDDIR}/polished.hq_isoforms.sorted.sam --dun-merge-5-shorter -o ${ALIGNEDDIR}/Collapsed/merged
   
## generate counts
get_abundance_post_collapse.py ${ALIGNEDDIR}/Collapsed/merged.collapsed ${PROCESSEDDIR}/Cluster/clustered.cluster_report.csv

## Filter away 5' degraded isoforms
filter_away_subset.py ${ALIGNEDDIR}/Collapsed/merged.collapsed

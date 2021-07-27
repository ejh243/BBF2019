#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p sq # submit to the serial queue
#SBATCH --time=72:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=4 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/MergeAndCluster.o
#SBATCH --error=LogFiles/MergeAndCluster.e
#SBATCH --job-name=MergeAndCluster

echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

module load Miniconda2
source activate isoseq
source ./Config/config.txt

echo "Changing Folder to: "
echo $DATADIR

cd $DATADIR/

## SMRT cells will be merged for subsquent steps 
cd $DATADIR/

## create file with list of all filenames
ls ${PROCESSEDDIR}/Refine/*.flnc.bam > ${PROCESSEDDIR}/Refine/flnc.fofn

## cluster data
mkdir -p ${PROCESSEDDIR}/Cluster 

## output will be split into 24 files so that polishing can be parallelized
isoseq3 cluster ${PROCESSEDDIR}/Refine/flnc.fofn ${PROCESSEDDIR}/Cluster/clustered.bam --verbose --use-qvs --log-file ${PROCESSEDDIR}/Cluster/merge.log --split-bam 24

## for polishing need to merge .xml files
XMLFILES=($(ls ${DATADIR}/*subreadset.xml))
dataset create --type SubreadSet ${DATADIR}/merged.subreadset.xml ${XMLFILES[@]}

#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=01:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mem=100G # 120 GB total 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/Cluster2_test-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/Cluster2_test-%j.err
#SBATCH --job-name=Cluster2_test

## bash script to automate merging and clustering of processed IsoSeq data
## Updated to use the latest Isoseq3 tools - version 4.3.0 (https://isoseq.how/getting-started.html)
## assumes isoseq3 has been installed/loaded 
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires .flnc.bam files are located in the ${MERGEDDIR} folder

# If multiple SMRTcells found in MERGEDDIR, then create a list of all <movie>.flnc.bam using
# This list can be used as input for cluster2 step

# Usage: sbatch Scripts/IsoSeqPipeline/mergeClusterIsoSeqSMRTcells.sh

source ./Config/config.txt
module load Miniconda3
source activate isoseq_tools   

echo "Looking for FLNC files in: ${MERGEDDIR}"
mkdir -p ${MERGEDDIR}/Cluster2


# Step 1: create fofn 
if [ ! -f ${MERGEDDIR}/flnc.fofn ]; then
    echo "Creating FLNC file list"
    find ${MERGEDDIR} -name "*.flnc.bam" > ${MERGEDDIR}/flnc.fofn
else
    echo "flnc.fofn exists - skipping"
fi
 

# Cluster2 - Cluster FLNC reads and generate transcripts
# Note - polish step not required in newer pipeline / cluster2 tool 

if [ ! -f ${MERGEDDIR}/Cluster2/clustered.bam ]; then
    echo "Running cluster2"

    isoseq cluster2 \
        ${MERGEDDIR}/flnc.fofn \
        ${MERGEDDIR}/Cluster2/clustered.bam \
        --singletons \
        --num-threads ${SLURM_CPUS_PER_TASK} \
        --log-file ${MERGEDDIR}/Cluster2/cluster2.log

else
    echo "Cluster2 output exists - skipping"
fi

echo "Done"


# End of script
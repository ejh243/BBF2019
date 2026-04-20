#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/Cluster2_test-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/Cluster2_test-%j.err
#SBATCH --job-name=Cluster2_test

## bash script to automate merging and clustering of processed IsoSeq data
## Updated to use the latest Isoseq3 tools - version 4.3.0 (https://isoseq.how/getting-started.html)
## assumes isoseq3 has been installed/loaded 
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires .flnc.bam files are located in the ${PROCESSEDDIR}/Refine folder

# If multiple SMRTcells found in PROCESSEDDIR, then create a list of all <movie>.flnc.bam using
# This list can be used as input for cluster2 step

# Usage: sbatch Scripts/IsoSeqPipeline/mergeClusterIsoSeqSMRTcells.sh

echo "Looking for FLNC files in: ${PROCESSEDDIR}"
mkdir -p ${PROCESSEDDIR}/Cluster2


# Step 1: create fofn 
if [ ! -f ${PROCESSEDDIR}/Refine/flnc.fofn ]; then
    echo "Creating FLNC file list"
    find ${PROCESSEDDIR}/Refine -name "*.flnc.bam" > ${PROCESSEDDIR}/Refine/flnc.fofn
else
    echo "flnc.fofn exists - skipping"
fi
 

# Cluster2 - Cluster FLNC reads and generate transcripts
# Note - polish step not required in newer pipeline / cluster2 tool 

if [ ! -f ${PROCESSEDDIR}/Cluster2/clustered_flnc_fofn.bam ]; then
    echo "Running cluster2"

    isoseq cluster2 \
        ${PROCESSEDDIR}/Refine/flnc.fofn \
        ${PROCESSEDDIR}/Cluster2/clustered_flnc_fofn.bam \
        --singletons \
        --log-file ${PROCESSEDDIR}/Cluster2/cluster2.log

else
    echo "Cluster output exists - skipping"
fi

echo "Done"


# End of script
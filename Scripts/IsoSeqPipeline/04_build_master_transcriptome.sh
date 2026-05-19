#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/BuildTranscriptome-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/BuildTranscriptome-%j.err
#SBATCH --job-name=BuildTranscriptome

## bash script to automate transcript clustering, read mapping and collapsing FLNC into unique isoforms 
## this script requires flnc.fofn are located in the ${PROCESSEDDIR}/FLNC folder

## do not store any sensitive data use config file to specify filepaths etc.

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/IsoSeqPipeline/04_build_master_transcriptome.sh

set -euo pipefail

echo "Starting Build Master Transcriptome job"
echo ""

## Load config and environment 
source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools 


## Set output dir
mkdir -p "${MERGEDDIR}/Clustered"
mkdir -p "${MERGEDDIR}/Aligned"
mkdir -p "${MASTERTRANSCRIPTOME}/Isoforms"


## Set relevant filepaths and output prefix 
fofn="${PROCESSEDDIR}/FLNC/flnc.fofn"

out_cluster="${MERGEDDIR}/Clustered/clustered"
out_mapped="${MERGEDDIR}/Aligned/mapped"
out_collapse="${MASTERTRANSCRIPTOME}/Isoforms/collapsed"


## Set safe temp directory
export TMPDIR="${MERGEDDIR}/tmp"
mkdir -p "${TMPDIR}"

echo "TMPDIR set to: $TMPDIR"
df -h "$TMPDIR"

## Step 1: IsoSeq Cluster2 - Cluster FLNC reads and generate transcripts
# if valid output file exists, skip step  
if [ -s "${out_cluster}.bam" ]; then 
    echo "IsoSeq Cluster2 output exists - skipping..."
else
    echo "Clustering similar transcripts across all SMRT cells..."

    # Check if required input file exists
    if [ ! -s "${fofn}" ]; then
        echo "ERROR: Missing or empty input file (flnc.fofn) in ${PROCESSEDDIR}/FLNC"
        exit 1
    fi

    # Run IsoSeq Cluster2
    isoseq cluster2 \
        "${fofn}" \
        "${out_cluster}.bam" \
        --singletons \
        --num-threads ${SLURM_CPUS_PER_TASK} \
        --log-file ${out_cluster}.log
fi


## Step 2: PacBio Minimap2 - Map all FLNC reads to reference genome 
# if valid output file exists, skip step  
if [ -s "${out_mapped}.bam" ]; then 
    echo "Aligned BAM exists - skipping..."
else
    echo "Aligning clustered transcripts against reference genome..."

    # Check if required input file exists
    if [ ! -s "${out_cluster}.bam" ]; then
        echo "ERROR: Missing or empty input file (clustered.bam)"
        exit 1
    fi

    # Run Minimap2 
    pbmm2 align \
    --preset ISOSEQ \
    --sort \
    --num-threads ${SLURM_CPUS_PER_TASK} \
    --log-file ${out_mapped}.log \
    ${REFGENOME} \
    ${out_cluster}.bam \
    ${out_mapped}.bam 
fi

## Step 3: IsoSeq Collapse - Collapse transcripts into unique isoforms 
# if valid output file exists, skip step  
if [ -s "${out_collapse}.gff" ]; then 
    echo "Collapsed isoforms exist - skipping..."
else
    echo "Collapsing mapped FLNC transcripts into unique isoforms..."

    # Check if required input file exists
    if [ ! -s "${out_mapped}.bam" ]; then
        echo "ERROR: Missing or empty input file (mapped.bam)"
        exit 1
    fi

    # Run IsoSeq Collapse
    isoseq collapse \
    --do-not-collapse-extra-5exons \
    --max-5p-diff 50 \
    --max-3p-diff 100 \
    --num-threads ${SLURM_CPUS_PER_TASK} \
    ${out_mapped}.bam \
    ${fofn} \
    ${out_collapse}.gff
fi

echo "Master Transcriptome complete"

rmdir -r "${TMPDIR}" # remove dir with temporary files 

## End of script
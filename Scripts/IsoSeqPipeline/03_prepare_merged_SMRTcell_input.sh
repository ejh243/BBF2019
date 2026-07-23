#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=01:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/PrepareMergedInput-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/PrepareMergedInput-%j.err
#SBATCH --job-name=PrepareMergedInput

## bash script to automate preparation of merging of individual processed SMRT cell data into file of filenames (fofn)
## this script requires processed flnc.bam files are located in the ${PROCESSEDDIR}/Refine folder
## SM header tags are updated to match SMRT cell IDs and a copy of the flnc.bams are placed in ${MERGEDDIR}
## flnc.fofn created for downstream analysis 

## do not store any sensitive data use config file to specify filepaths etc.

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/IsoSeqPipeline/03_prepare_merged_SMRTcell_input.sh

set -euo pipefail

echo "Preparing individual SMRT cell FLNC data as required for subsequent merged analysis..."

## Load config and environment 
source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools 


# Set output dir
mkdir -p "${MERGEDDIR}/PreparedFLNC"


## Check how many FLNC BAM files exist
echo "Scanning for FLNC BAM files in: ${PROCESSEDDIR}/Refine..."

bam_files=("${PROCESSEDDIR}"/Refine/*.flnc.bam)

if [ ${#bam_files[@]} -eq 0 ]; then
    echo "ERROR: No FLNC BAM files found in ${PROCESSEDDIR}/Refine"
    exit 1
fi

echo "Found ${#bam_files[@]} SMRT cell(s)"


## Prepare all SMRT cell runs 
for bam in "${bam_files[@]}"; do
    run=$(basename "$bam" .flnc.bam)
    renamed_out="${MERGEDDIR}/PreparedFLNC/${run}.flnc.bam"

    echo "Processing: $run"

    # Skip if already done
    if [ -f "$renamed_out" ]; then
        echo "Output already exists - skipping"
        continue
    fi

    # Replace SM tag in BAM header
    echo "Renaming SM tag in BAM header..."

    samtools view -H "$bam" \
    | sed -E "/^@RG/ s/SM:[^ \t]+/SM:SMRT_${run}/" \
    | tee /dev/stderr \
    | samtools reheader - "$bam" > "$renamed_out"

    # Sanity check
    if ! samtools view -H "$renamed_out" | grep -q "SM:SMRT_${run}"; then
        echo "ERROR: SM tag update failed for $run"
        exit 1
    fi


    # Generate PBI index file
    echo "Generating PBI index..."
    pbindex "$renamed_out"
    echo 

done

# [optional] can run quick check using: samtools view -H "$out" | grep '^@RG'


# Create file of filenames  
echo  
if [ ! -f ${MERGEDDIR}/PreparedFLNC/flnc.fofn ]; then
    echo "Creating file of filenames (fofn)..."
    find ${MERGEDDIR}/PreparedFLNC/ -name "*.flnc.bam" > ${MERGEDDIR}/PreparedFLNC/flnc.fofn
else
    echo "flnc.fofn exists - skipping"
fi

echo 
echo "Individual SMRT cell data (FLNC reads) prepared. Ready for merged analysis."

# End of script 
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
#SBATCH --output=/lustre/home/vs455/LogFiles/RenameSample-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/RenameSample-%j.err
#SBATCH --job-name=RenameSample

## bash script to automate renaming of samples (SM tag) in processed BAM files
## do not store any sensitive data use config file to specify filepaths etc.


# Usage: sbatch Scripts/IsoSeqPipeline/changeBAMHeaders.sh

source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools  


# Set output dir
mkdir -p "${MERGEDDIR}"

# Automate to all SMRT cell runs 
for bam in "${PROCESSEDDIR}"/Refine/*.flnc.bam; do
    run=$(basename "$bam" .flnc.bam)
    renamed_out="${MERGEDDIR}/${run}.renamed_flnc.bam"

    echo "Processing $run"

    # Skip if already done
    if [ ! -f "$renamed_out" ]; then
        echo "Renaming SM header..."

        # Replace SM tag in BAM header
        samtools view -H "$bam" \
        | sed -E "s/SM:[^[:space:]]+/SM:SMRT_${run}/" \
        | samtools reheader - "$bam" > "$renamed_out"

        # Generate PBI index file
        echo "Generating PBI index..."
        pbindex "$renamed_out"

    else
        echo "Renamed BAM file exists - skipping"
    fi
done


# Quick check using: 
# samtools view -H "$out" | grep '^@RG'


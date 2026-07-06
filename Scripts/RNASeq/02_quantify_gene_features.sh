#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=01:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=4 
#SBATCH --mem=8G 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lfs1i3/projects/e6e/LogFiles/FeatureCounts_ExSR-%j.out 
#SBATCH --error=/lfs1i3/projects/e6e/LogFiles/FeatureCounts_ExSR-%j.err 
#SBATCH --job-name=FeatureCounts


## bash script to quantify gene features based on read count
## this script requires aligned BAM files are located in the ${RNASEQDIR}/aligned_reads folder

## do not store any sensitive data use config file to specify filepaths etc.

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/RNASeq/02_quantify_gene_features.sh

set -euo pipefail

echo "Starting Gene-level quantification job..."
echo


## Load required software and configurations 
source ./Config/config_v2.txt

source ~/miniconda3/etc/profile.d/conda.sh  # in place of module load Miniconda3
conda activate rnaseq_tools


## Set variables  
THREADS=${SLURM_CPUS_PER_TASK:-16}

ALIGNEDDIR="${RNASEQDIR}/aligned_reads_primary"
OUTPUT="${RNASEQDIR}/gene_counts_primary_bam_s2.txt"


## Check input BAM files exist
num_bams=$(ls "${ALIGNEDDIR}"/*.bam 2>/dev/null | wc -l)

if [[ $num_bams -eq 0 ]]; then
    echo "ERROR: No BAM files found in ${ALIGNEDDIR}"
    exit 1
fi

echo "Found $num_bams BAM files in ${ALIGNEDDIR}"
echo


## Run featureCounts 
if [[ -f "$OUTPUT" ]]; then
    echo "Gene count file already exists, skipping featureCounts..."
else
    echo "Running featureCounts..."

    # Using default settings for paired short reads
    featureCounts \
      -T "$THREADS" \
      -p \
      -s 2 \
      -t exon \
      -g gene_id \
      -a "${GENCODEGTF}" \
      -o "$OUTPUT" \
      "${ALIGNEDDIR}"/*.bam 
      
    echo "featureCounts completed"
fi

echo
echo "Job completed..."

# End of script
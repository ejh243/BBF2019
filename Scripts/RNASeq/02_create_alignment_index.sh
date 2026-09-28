#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=01:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation
#SBATCH --mem=64G # Memory usage 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lfs1i3/projects/e6e/LogFiles/Create_Index-%j.out
#SBATCH --error=/lfs1i3/projects/e6e/LogFiles/Create_Index-%j.err
#SBATCH --job-name=Create_Index


## script to create STAR Index and Salmon Index ahead of alignment step
## STAR Index requires ${REFGENOME} and annotation ${GTF} to be located in ${RESOURCESDIR}
## Salmon Index requires ${TRANSCRIPTSFA} to be located in ${RESOURCESDIR}

# Usage: sbatch Scripts/RNASeq/02_create_alignment_index.sh

set -euo pipefail

echo "Starting STAR and Salmon Index creation job..."
echo

## Load required software and configurations 
source ./Config/config_v2.txt

source ~/miniconda3/etc/profile.d/conda.sh  # in place of module load Miniconda3
conda activate rnaseq_tools

# Set output dir
STARDIR="${RESOURCESDIR}/STARIndex"
SALMONDIR="${RESOURCESDIR}/SalmonIndex"
mkdir -p "$STARDIR" 
mkdir -p "$SALMONDIR" 


## Run STAR on genomeGenerate mode
if [[ -f "${STARDIR}/sjdbInfo.txt" ]]; then
    echo "STAR Index already exists. Skipping index creation."
else
    echo "Building STAR index..."

    STAR \
      --runThreadN ${SLURM_CPUS_PER_TASK} \
      --runMode genomeGenerate \
      --genomeDir ${STARDIR} \
      --genomeFastaFiles ${REFGENOME} \
      --sjdbGTFfile ${GTF} \
      --sjdbOverhang 99

    echo "Output written to:"
    echo "$STARDIR"
fi

echo 

## Create Salmon index for custom transcriptome 
if [[ -f "${SALMONDIR}/info.json" ]]; then
    echo "Salmon index already exists. Skipping index creation."
else
    echo "Building Salmon index..."

    salmon index \
        -t "${TRANSCRIPTSFA}" \
        -i "$SALMONDIR" \
        -p ${SLURM_CPUS_PER_TASK}

    echo "Output written to:"
    echo "$SALMONDIR"
fi

echo 
echo "STAR and Salmon Indexes creation job completed."

# End of script
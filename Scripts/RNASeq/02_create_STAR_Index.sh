#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=05:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation
#SBATCH --mem=64G # Memory usage 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lfs1i3/projects/e6e/LogFiles/STAR_Index-%j.out
#SBATCH --error=/lfs1i3/projects/e6e/LogFiles/STAR_Index-%j.err
#SBATCH --job-name=STAR_Index


## script to create STAR index ahead of alignment 

# Usage: sbatch Scripts/RNASeq/create_STAR_Index.sh

set -euo pipefail


## Load required software and configurations 
source ./Config/config_v2.txt

source ~/miniconda3/etc/profile.d/conda.sh  # in place of module load Miniconda3
conda activate rnaseq_tools


# Set output dir
STARDIR="${RESOURCESDIR}/STARIndex"
mkdir -p "$STARDIR"


## Run STAR on genomeGenerate mode
STAR \
  --runThreadN ${SLURM_CPUS_PER_TASK} \
  --runMode genomeGenerate \
  --genomeDir ${STARDIR} \
  --genomeFastaFiles ${REFGENOME} \
  --sjdbGTFfile ${GENCODEGTF} \
  --sjdbOverhang 99

# End of script
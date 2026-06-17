#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=05:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/STAR_Index-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/STAR_Index-%j.err
#SBATCH --job-name=STAR_Index

## script to create STAR index ahead of alignment 

# Usage: sbatch Scripts/RNASeq/create_STAR_Index.sh

set -euo pipefail


## Load required software and configurations 
source ./Config/config.txt

module load STAR
module load Miniconda3
source activate rnaseq_tools


# Set output dir
STARDIR="${RESOURCESDIR}/STARIndex"
mkdir -p "$STARDIR"


## Run STAR on genomeGenerate mode
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir ${STARDIR} \
  --genomeFastaFiles ${REFGENOME} \
  --sjdbGTFfile ${GENCODEGTF} \
  --sjdbOverhang 99

# End of script
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=05:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation
#SBATCH --mem=64G # Memory usage 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lfs1i3/projects/e6e/LogFiles/STAR-%A_%a.out 
#SBATCH --error=/lfs1i3/projects/e6e/LogFiles/STAR-%A_%a.out 
#SBATCH --job-name=STAR_Alignment
#SBATCH --array=0-47%3 ## runs multiple jobs with 3 at any one time 

## bash script to automate preprocessing of paired short read data 
## Parallelisation: Uses a SLURM job array (one SMRT cell per task)
    ## Configure via --array=0-N%M where N = samples-1 and M = max concurrent jobs

## script to automate gene-level alignment (using STAR) 

## do not store any sensitive data use config file to specify filepaths etc. 

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/RNASeq/03_perform_read_alignment.sh 

set -euo pipefail

echo "Starting Read Alignment job..."
echo


## Load required software and configurations 
source ./Config/config_v2.txt

source ~/miniconda3/etc/profile.d/conda.sh  
conda activate rnaseq_tools

## Relevant directories
TRIMDIR="${RNASEQDIR}/reads_trimmed"
INDEX="${RESOURCESDIR}/STARIndex"
ALIGNEDRNA="${RNASEQDIR}/aligned_reads"

## Extract sampleName from TRIMDIR 
mapfile -t SAMPLES < <(
    find "${TRIMDIR}" -maxdepth 1 -name "*val_1.f*q.gz" |
    xargs -n1 basename |
    sed -E 's/_R1_.*_val_1\.f.*gz$//' |
    sort -u
)

sampleName="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

## Safety checks 
if [[ ${SLURM_ARRAY_TASK_ID} -ge ${#SAMPLES[@]} ]]; then
    echo "ERROR: Array index ${SLURM_ARRAY_TASK_ID} exceeds number of samples (${#SAMPLES[@]})"
    exit 1
fi


## Run STAR alignment 
echo
echo "Preparing STAR alignment..."

## Locate trimmed reads dynamically
star_f1=$(ls "${TRIMDIR}/${sampleName}"*val_1.f*q.gz 2>/dev/null | head -n 1)
star_f2=$(ls "${TRIMDIR}/${sampleName}"*val_2.f*q.gz 2>/dev/null | head -n 1)


## Validate
[[ -f "$star_f1" ]] || { echo "ERROR: Trimmed R1 not found"; exit 1; }
[[ -f "$star_f2" ]] || { echo "ERROR: Trimmed R2 not found"; exit 1; }

echo "Using trimmed reads:"
echo "$star_f1"
echo "$star_f2"


## Align with STAR using star index created previously 
if [[ ! -f "${ALIGNEDRNA}/${sampleName}Log.final.out" ]]; then
    echo "Starting STAR alignment..."

    STAR \
        --genomeDir "${INDEX}" \
        --runThreadN "$THREADS" \
        --readFilesIn "$star_f1" "$star_f2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${ALIGNEDRNA}/${sampleName}" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard

    echo "STAR alignment complete."
else
    echo "STAR already completed for ${sampleName}, skipping"
fi

echo
echo "Pipeline completed for ${sampleName}"
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=01:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --cpus-per-task=16 # specify number of threads per task 
#SBATCH --mem=64G # Memory usage 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # enter email address
#SBATCH --output=/lfs1i3/projects/e6e/LogFiles/Salmon-%A_%a.out 
#SBATCH --error=/lfs1i3/projects/e6e/LogFiles/Salmon-%A_%a.err 
#SBATCH --job-name=SalmonQuant
#SBATCH --array=0-47%3 ## runs multiple jobs with 3 at any one time 

## bash script to quantify transcript isoforms using salmon 
## do not store any sensitive data use config file to specify filepaths etc.

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/RNASeq/05_quantify_transcripts.sh

set -euo pipefail

echo "Starting Transcript quantification job..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo 


## Load required software and configurations 
source ./Config/config_v2.txt

source ~/miniconda3/etc/profile.d/conda.sh  # in place of module load Miniconda3
conda activate rnaseq_tools


## Set variables  
THREADS=${SLURM_CPUS_PER_TASK:-16}

TRIMDIR="${RNASEQDIR}/reads_trimmed"
INDEX="${RESOURCESDIR}/SalmonIndex"
OUTDIR="${RNASEQDIR}/transcript_quant"

mkdir -p "$OUTDIR" 


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
read1=$(ls "${TRIMDIR}/${sampleName}"*val_1.f*q.gz 2>/dev/null | head -n 1)
read2=$(ls "${TRIMDIR}/${sampleName}"*val_2.f*q.gz 2>/dev/null | head -n 1)


## Validate
[[ -f "$read1" ]] || { echo "ERROR: Trimmed R1 not found"; exit 1; }
[[ -f "$read2" ]] || { echo "ERROR: Trimmed R2 not found"; exit 1; }

echo "Using trimmed reads:"
echo "$read1"
echo "$read2"

mkdir -p "${OUTDIR}/${sampleName}"

## Run Salmon transcript quantification
if [[ -f "${OUTDIR}/${sampleName}/quant.sf" ]]; then
    echo "Salmon Quantification file already exists for ${sampleName}, skipping..."
else
    echo "Running Salmon Quantification..."

    # Using default settings for paired short reads
    salmon quant \
    -i "${INDEX}" \
    -l A \
    -1 "${read1}" \
    -2 "${read2}" \
    -p "$THREADS" \
    -o "${OUTDIR}/${sampleName}"

    echo "Salmon quantification complete for ${sampleName}"
    echo "Output written to:"
    echo "${OUTDIR}"
fi

echo
echo "Job completed..."

# End of script




#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --cpus-per-task=16 # specify number of threads per task 
#SBATCH --mem=64G # Memory usage 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # enter email address
#SBATCH --output=/lfs1i3/projects/e6e/LogFiles/Salmon_ExSR-%A_%a.out 
#SBATCH --error=/lfs1i3/projects/e6e/LogFiles/Salmon_ExSR-%A_%a.err 
#SBATCH --job-name=SalmonQuant
#SBATCH --array=0-1%3 ## runs multiple jobs with 3 at any one time 

## bash script to quantify transcript isoforms using salmon 
## do not store any sensitive data use config file to specify filepaths etc.

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/RNASeq/04_quantify_transcripts.sh

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

FASTQDIR="${RNASEQDIR}/reads_trimmed"
INDEX="${RESOURCESDIR}/SalmonIndex"
OUTPUTDIR="${RNASEQDIR}/salmonQuant"

mkdir -p "$OUTPUTDIR" 


## Check input FASTQ files exist
num_fastqs=$(ls "${FASTQDIR}"/*R1*.fq.gz 2>/dev/null | wc -l)

if [[ $num_fastqs -eq 0 ]]; then
    echo "ERROR: No trimmed FASTQ files found in ${FASTQDIR}"
    exit 1
fi
echo "Found $num_fastqs paired-end FASTQ files."
echo


## Assign sample per array task 
FASTQS=("${FASTQDIR}"/*_val_1.f*q.gz)

R1="${FASTQS[$SLURM_ARRAY_TASK_ID]}"
sampleName= # Strip everything from R1 onward

R2=""
10_TAGCTT_L007_R1_001_val_1.fq.gz
10_TAGCTT_L007_R2_001_val_2.fq.gz


OUT="${OUTPUTDIR}/${sampleName}"

echo "Processing sample: ${sampleName}"
echo "R1: $R1"
echo "R2: $R2"
echo


## Run Salmon transcript quantification
if [[ -f "${OUT}/quant.sf" ]]; then
    echo "Salmon Quantification file already exists for ${SAMPLE}, skipping..."
else
    echo "Running Salmon Quantification..."

    # Using default settings for paired short reads
    salmon quant \
    -i "$INDEX" \
    -l A \
    -1 "$R1" \
    -2 "$R2" \
    -p "$THREADS" \
    -o "$OUT"

    echo "Salmon quantification complete."
    echo "Output written to:"
    echo "$OUT"
fi


echo
echo "Job completed..."

# End of script




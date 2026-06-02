#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # enter email address
#SBATCH --output=/lustre/home/vs455/LogFiles/PreprocessEpiGABA-%A_%a.out 
#SBATCH --error=/lustre/home/vs455/LogFiles/PreprocessEpiGABA-%A_%a.err 
#SBATCH --job-name=PreprocessEpiGABA
#SBATCH --array=0-30%6 ## runs multiple jobs with 6 at any one time 

## bash script to automate preprocessing of paired short read data 
## Parallelisation: Uses a SLURM job array (one SMRT cell per task)
    ## Configure via --array=0-N%M where N = samples-1 and M = max concurrent jobs

## do not store any sensitive data use config file to specify filepaths etc. 

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/RNASeq/01_preprocess_ShortReads.sh

set -euo pipefail

echo "Starting RNA-Seq preprocessing job"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo ""


## Load required software and configurations 
source ./Config/config.txt

module load FastQC
module load Miniconda3
source activate rnaseq_tools

# output software versions 
echo "software tools used"
trim_galore --version
fastqc --version


## Output directories 
TRIMDIR="${RNASEQDIR}/reads_trimmed"
FASTQC_RAW="${RNASEQDIR}/fastqc_raw"
FASTQC_TRIMMED="${RNASEQDIR}/fastqc_trimmed"

mkdir -p "${TRIMDIR}"
mkdir -p "${FASTQC_RAW}"
mkdir -p "${FASTQC_TRIMMED}"


## Locate all input FASTQ files 
mapfile -t ALL_FASTQS < <(
    find "${SHORTREADS}" \
        -maxdepth 1 \
        -name "*.fastq.gz" \
        | sort
)

echo "Total FASTQ.GZ files found: ${#ALL_FASTQS[@]}"

FQFILES=($(find ${SHORTREADS} -maxdepth 1 -name '*.fastq.gz' ))

# Build sample list using R1 files only
mapfile -t FQFILES < <(
    find "${SHORTREADS}" \
        -maxdepth 1 \
        -name "*.R1.fastq.gz" \
        | sort
)

echo "Number of samples found: ${#FQFILES[@]}"
echo

# Safety check 
if [[ ${SLURM_ARRAY_TASK_ID} -ge ${#FQFILES[@]} ]]; then
    echo "ERROR: Array index ${SLURM_ARRAY_TASK_ID} exceeds number of samples (${#FQFILES[@]})."
    exit 1
fi

# Select sample for this array
f1="${FQFILES[$SLURM_ARRAY_TASK_ID]}"
f2="${f1/.R1.fastq.gz/.R2.fastq.gz}"

sampleName=$(basename "${f1}" .R1.fastq.gz)

echo "Processing sample: ${sampleName}"
echo "R1 FASTQ: ${f1}"
echo "R2 FASTQ: ${f2}"
echo


# Verify paired FASTQ exists 
if [[ ! -f "${f2}" ]]; then
    echo "ERROR: Mate pair not found."
    echo "Expected file:"
    echo "${f2}"
    exit 1
fi

echo "Paired FASTQ files verified."
echo


# Step 1: Run FASTQC on Raw Reads
echo "Running FastQC on raw reads..."

fastqc "${f1}" "${f2}" \
    --threads 16 \
    --outdir "${FASTQC_RAW}"


# Step 2: Run Trim Galore 
echo "Running Trim Galore..."

trimmed_f1="${TRIMDIR}/${sampleName}.R1.fastq.gz"
trimmed_f2="${TRIMDIR}/${sampleName}.R2.fastq.gz"

trim_galore --paired ${f1} ${f2} --fastqc -o ${TRIMDIR}

if [[ -s "${trimmed_f1}" ]] || [[-s "${trimmed_f2}" ]]; then
echo "ERROR: Trimmed FASTQ files not found"
exit 1
fi


# Step 3: Move FASTQC results into relevant folder 
echo "Running FastQC on trimmed reads..."

fastqc "${trimmed_f1}" "${trimmed_f2}" \
    --threads 16 \
    --outdir "${FASTQC_TRIMMED}"

echo
echo "Completed preprocessing for ${sampleName}"

# End of script 
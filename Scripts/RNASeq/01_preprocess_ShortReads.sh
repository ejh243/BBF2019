#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # enter email address
#SBATCH --output=/lustre/home/vs455/LogFiles/PreprocessShortReads-%A_%a.out 
#SBATCH --error=/lustre/home/vs455/LogFiles/PreprocessShortReads-%A_%a.err 
#SBATCH --job-name=PreprocessShortReads
#SBATCH --array=0-19%5 ## runs multiple jobs with 5 at any one time 

## bash script to automate preprocessing of paired short read data 
## Parallelisation: Uses a SLURM job array (one SMRT cell per task)
    ## Configure via --array=0-N%M where N = samples-1 and M = max concurrent jobs

## do not store any sensitive data use config file to specify filepaths etc. 
## please provide either $SHORTREADS or $SHORTREAD_LIST in the config file
#   Assumes paired-end FASTQs are named: sample.R1.fastq.gz & sample.R2.fastq.gz
#   If a different naming convention is used, modify lines flagged with "# must match filename pattern" 

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

# Validate input configuration 
if [[ -n "${SHORTREADS:-}" && -n "${SHORTREAD_LIST:-}" ]]; then
    echo "ERROR: Specify either SHORTREADS or SHORTREAD_LIST in config file, not both"
    exit 1
fi

if [[ -z "${SHORTREADS:-}" && -z "${SHORTREAD_LIST:-}" ]]; then
    echo "ERROR: Specify either SHORTREADS or SHORTREAD_LIST in config file"
    exit 1
fi


## Read FASTQ files
if [[ -n "${SHORTREAD_LIST:-}" ]]; then # If list is provided, check if files exist
    if [[ ! -f "${SHORTREAD_LIST}" ]]; then
        echo "ERROR: FASTQ list file not found:"
        echo "${SHORTREAD_LIST}"
        exit 1
    fi

    echo "Reading FASTQ files from list:"
    echo "${SHORTREAD_LIST}"
    mapfile -t ALL_FASTQS < "${SHORTREAD_LIST}"

else # If filepath provided, read all FASTQ files in provided directory 
    echo "Searching FASTQ files in:"
    echo "${SHORTREADS}"
    mapfile -t ALL_FASTQS < <(find "${SHORTREADS}" -maxdepth 1 -name "*.fastq.gz" | sort)
fi

echo "Total FASTQ files found: ${#ALL_FASTQS[@]}"
echo


## Build sample list using R1 FASTQ files
mapfile -t FQFILES < <(
    printf "%s\n" "${ALL_FASTQS[@]}" |
    grep '_R1_001\.fastq\.gz$' |   # must match filename pattern 
    sort
)

echo "Number of samples found: ${#FQFILES[@]}"
echo


# Safety checks 
if [[ ${#FQFILES[@]} -eq 0 ]]; then
    echo "ERROR: No R1 FASTQ files found."
    exit 1
fi

if [[ ${SLURM_ARRAY_TASK_ID} -ge ${#FQFILES[@]} ]]; then
    echo "ERROR: Array index ${SLURM_ARRAY_TASK_ID} exceeds number of samples (${#FQFILES[@]})."
    exit 1
fi


# Select sample for this array
f1="${FQFILES[$SLURM_ARRAY_TASK_ID]}"
f2="${f1/_R1_001.fastq.gz/_R2_001.fastq.gz}" # must match filename pattern 

sampleName=$(basename "${f1}" _R1_001.fastq.gz) # must match filename pattern 

echo "Processing sample: ${sampleName}"
echo "R1 FASTQ: ${f1}"
echo "R2 FASTQ: ${f2}"
echo

## Verify paired FASTQ exists 
if [[ ! -f "${f1}" ]]; then
    echo "ERROR: R1 FASTQ file not found."
    echo "Expected file:"
    echo "${f1}"
    exit 1
fi

if [[ ! -f "${f2}" ]]; then
    echo "ERROR: R2 FASTQ file not found."
    echo "Expected file:"
    echo "${f2}"
    exit 1
fi

echo "Paired FASTQ files verified."
echo


# Step 1: Run FASTQC on Raw Reads
echo "Running FastQC on raw reads..."
#fastqc "${f1}" "${f2}"  --threads 16  --outdir "${FASTQC_RAW}"


# Step 2: Run Trim Galore 
echo "Running Trim Galore..."
echo "followed by FastQC on trimmed reads..."
trim_galore --paired  "${f1}" "${f2}" --fastqc -o "${TRIMDIR}"


# Step 3: Move FASTQC reports to relevant folder 
echo "Moving FastQC reports..."
mv "${TRIMDIR}"/*fastqc.* "${FASTQC_TRIMMED}/" 2>/dev/null || true

echo
echo "Completed preprocessing for ${sampleName}"

# End of script 
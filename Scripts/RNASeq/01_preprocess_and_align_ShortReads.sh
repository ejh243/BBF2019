#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --cpus-per-task=16 # specify number of threads per task 
#SBATCH --mem=64G # Memory usage 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # enter email address
#SBATCH --output=/lfs1i3/projects/e6e/LogFiles/AlignShortReads-%A_%a.out 
#SBATCH --error=/lfs1i3/projects/e6e/LogFiles/AlignShortReads-%A_%a.err 
#SBATCH --job-name=AlignShortReads
#SBATCH --array=0-47%3 ## runs multiple jobs with 3 at any one time 

## bash script to automate preprocessing of paired short read data 
## Parallelisation: Uses a SLURM job array (one SMRT cell per task)
    ## Configure via --array=0-N%M where N = samples-1 and M = max concurrent jobs

## do not store any sensitive data use config file to specify filepaths etc. 

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/RNASeq/01_preprocess_and_align_ShortReads.sh

set -euo pipefail

echo "Starting RNA-Seq preprocessing job"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo ""

## Load required software and configurations 
source ./Config/config_v2.txt

source ~/miniconda3/etc/profile.d/conda.sh  # in place of module load Miniconda3
conda activate rnaseq_tools

# output software versions 
echo "software tools used"
trim_galore --version
fastqc --version
STAR --version 

THREADS=${SLURM_CPUS_PER_TASK:-16}


## Output directories 
TRIMDIR="${RNASEQDIR}/reads_trimmed"
FASTQC_RAW="${RNASEQDIR}/fastqc_raw"
FASTQC_TRIMMED="${RNASEQDIR}/fastqc_trimmed"
ALIGNEDRNA="${RNASEQDIR}/aligned_reads"

mkdir -p "$TRIMDIR" "$FASTQC_RAW" "$FASTQC_TRIMMED" "$ALIGNEDRNA"


## Locate all input FASTQ files 
echo "Searching FASTQ files in: ${SHORTREADS}"

mapfile -t ALL_FASTQS < <(find "${SHORTREADS}" -maxdepth 1 -name "*.fastq.gz" | sort)

echo "Total FASTQ files found: ${#ALL_FASTQS[@]}"
echo


## Identify damples (R1 filename based detection)
mapfile -t FQFILES < <(
    printf "%s\n" "${ALL_FASTQS[@]}" |
    grep -Ei 'r1.*\.fastq\.gz$' | 
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


## Get sample specific files
f1="${FQFILES[$SLURM_ARRAY_TASK_ID]}" 
f2=$(echo "$f1" | sed -E 's/R1/R2/; s/r1/r2/')  # R1 -> R2 substitution (case insensitive)
sampleName=$(basename "$f1" | sed -E 's/[._-]?[Rr]1.*\.fastq\.gz//') # Strip everything from R1 onward


echo "Processing sample: ${sampleName}"
echo "R1 FASTQ: ${f1}"
echo "R2 FASTQ: ${f2}"
echo


## Verify paired FASTQ exists 
if [[ ! -f "${f1}" ]]; then
    echo "ERROR: R1 FASTQ file not found."
    echo "Expected file:${f1}"
    exit 1
fi

if [[ ! -f "${f2}" ]]; then
    echo "ERROR: R2 FASTQ file not found."
    echo "Expected file:${f2}"
    exit 1
fi

echo "Paired FASTQ files verified."
echo


# Step 1: Run FASTQC on Raw Reads
if ! ls "${FASTQC_RAW}/${sampleName}"*fastqc.zip > /dev/null 2>&1; then
    echo "Running FastQC on raw reads..."
    fastqc "$f1" "$f2" --threads $THREADS --outdir "$FASTQC_RAW"
else
    echo "FastQC raw already exists, skipping..."
fi


# Step 2: Run Trim Galore 
echo
if ! ls "${TRIMDIR}/${sampleName}"*val_1.f*q.gz >/dev/null 2>&1; then
    echo "Running Trim Galore, followed by FastQC... "
    trim_galore --paired "$f1" "$f2" --fastqc -o "$TRIMDIR"
    mv "${TRIMDIR}"/*fastqc.* "${FASTQC_TRIMMED}/" 2>/dev/null || true

else
    echo "Trimmed files already exist, skipping..."
fi

echo "Completed preprocessing for ${sampleName}"



# Step 3: Run STAR alignment 
echo
echo "Preparing STAR alignment..."

# Locate trimmed reads dynamically
star_f1=$(ls "${TRIMDIR}/${sampleName}"*val_1.f*q.gz 2>/dev/null | head -n 1)
star_f2=$(ls "${TRIMDIR}/${sampleName}"*val_2.f*q.gz 2>/dev/null | head -n 1)

# Validate
[[ -f "$star_f1" ]] || { echo "ERROR: Trimmed R1 not found"; exit 1; }
[[ -f "$star_f2" ]] || { echo "ERROR: Trimmed R2 not found"; exit 1; }

echo "Using trimmed reads:"
echo "$star_f1"
echo "$star_f2"


## align with STAR using GENCODE v48 star index
if [[ ! -f "${ALIGNEDRNA}/${sampleName}Log.final.out" ]]; then
    echo "Starting STAR alignment..."

    STAR \
        --genomeDir "${RESOURCESDIR}/STARIndex" \
        --runThreadN $THREADS \
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

# End of script 
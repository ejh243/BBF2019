#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # enter email address
#SBATCH --output=/lustre/home/vs455/LogFiles/PreprocessIsoseq-%A_%a.out 
#SBATCH --error=/lustre/home/vs455/LogFiles/PreprocessIsoseq-%A_%a.err 
#SBATCH --job-name=PreprocessIsoseq
#SBATCH --array=0-34%7 ## runs multiple jobs with 7 at any one time 

## bash script to automate preprocessing of individual SMRT cells using PacBio IsoSeq tools
## Parallelisation: Uses a SLURM job array (one SMRT cell per task)
    ## Configure via --array=0-N%M where N = samples-1 and M = max concurrent jobs
## this script requires .subread.bam, .subreads.bam.pbi, and .subreadset.xml files are located in the DATADIR

## do not store any sensitive data use config file to specify filepaths etc. 

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/IsoSeqPipeline/01_preprocess_IsoSeqSMRTcells.sh

set -euo pipefail

echo "Starting IsoSeq preprocessing job"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"

source ./Config/config.txt

echo "Changing Folder to DATADIR: "
echo "${DATADIR}"
cd ${DATADIR}


## this command can be used to process all relevant files in DATADIR 
## if only a subset need to be processed provide list in FilesToProcess.txt and hash out line below.
samples=( *.subreads.bam )

echo "Number of samples found: " ${#samples[@]}

## run initial steps on each SMRT cell individually
sample=${DATADIR}/${samples[${SLURM_ARRAY_TASK_ID}]}

echo "Sample selected: ${sample}"

if [ ! -f "${sample}" ]; then
    echo "Input file not found, exiting:"
    echo "${sample}"
    exit 1
fi

basename=$(basename "${sample}" .subreads.bam)

## create output directories if required 
mkdir -p ${PROCESSEDDIR}/CCS
mkdir -p ${PROCESSEDDIR}/Lima
mkdir -p ${PROCESSEDDIR}/Refine

echo "Changing folder to Scripts directory: "
echo "${SCRIPTSDIR}/IsoSeqPipeline"
cd ${SCRIPTSDIR}/IsoSeqPipeline

## Load required software 
module load Miniconda3
source activate isoseq_tools   

# output software versions 
echo "software tools used"
isoseq --version
ccs --version
lima --version


## Assign variables for input and output file prefix 
ccs_output="${PROCESSEDDIR}/CCS/${basename}.ccs"
lima_output="${PROCESSEDDIR}/Lima/${basename}.fl"
refine_output="${PROCESSEDDIR}/Refine/${basename}.flnc"


## Main IsoSeq processing script ## 
THREADS=16

echo "Processing sample:"
echo "${basename}"

# Step 1: CCS - Generate circular consensus sequences (ccs) from subreads
# default min-passes is 3 full length subreads 
if [ -s "${ccs_output}.bam" ]; then # if valid output file exists, skip step  
    echo "CCS output file exists - skipping..."
else
    echo "Running Circular Consensus Sequence calling"

    # Check if required input file exists
    if [ ! -s "${sample}" ]; then
        echo "ERROR: Missing or empty input file for CCS (Raw data)"
        exit 1
    fi

    # Run CCS 
    ccs "${sample}" \
        "${ccs_output}.bam" \
        --min-rq 0.9 \
        --min-passes 1 \
        --num-threads $THREADS \
        --report-file ${ccs_output}_report.txt
fi


# Step 2: Lima - Remove cDNA primers and demultiplexing barcoded data 
# Lima appends primer-specific suffixes to ${lima_output}.<primer_5p--primer_3p>.bam
# Hence * here is used as a placeholder
if ls ${lima_output}.*.bam 1> /dev/null 2>&1; then # if valid output file exists, skip step 
    echo "Lima output file exists - skipping..."
else
    echo "Running Primer removal and Demultiplexing"

    # Check if required input file exists
    if [ ! -s "${ccs_output}.bam" ]; then
        echo "ERROR: Missing or empty input file for Lima (CCS output)"
        exit 1
    fi

    # Run Lima 
    lima --isoseq \
        --peek-guess --dump-clips \
        --num-threads $THREADS \
        ${ccs_output}.bam \
        ${PRIMERSEQ} \
        ${lima_output}.bam
fi


# Step 3: Refine - Remove polyA and concatemers from FL reads and generate FLNC transcripts
if [ -s "${refine_output}.bam" ]; then # if valid output file exists, skip step
    echo "Refine output file exists - skipping..."
else
    echo "Running Isoseq Refine"

    # Check if required input file exists
    if ! ls ${lima_output}.*.bam 1> /dev/null 2>&1; then
        echo "ERROR: Missing or empty input file for Refine (Lima output)"
        exit 1
    fi

    # Run IsoSeq Refine  
    isoseq refine \
        --require-polya \
        ${lima_output}.*.bam \
        ${PRIMERSEQ} \
        ${refine_output}.bam

    # Final output validation
    if [ ! -s "${refine_output}.bam" ]; then
        echo "ERROR: Refine failed — output not created"
        exit 1
    fi
fi


echo "Processing completed for:"
echo "${basename}"
# End of script 

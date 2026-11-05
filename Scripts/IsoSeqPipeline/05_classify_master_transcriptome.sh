#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=06:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/ClassifyIsoforms-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/ClassifyIsoforms-%j.err
#SBATCH --job-name=ClassifyIsoforms

## bash script to automate classification of collapsed isoforms using SQANTI

## script requires reference annotation GTF, collapsed isoforms GFF, FLNC count TXT, reference genome FASTA 
## and short read junctions coverage from STAR alignment to same reference (SJ.out.tab)

## do not store any sensitive data use config file to specify filepaths etc.

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/IsoSeqPipeline/05_classify_master_transcriptome.sh

set -euo pipefail

echo "Starting SQANTI3 Isoform classification job..."
echo 

## Load config and environment 
source ./Config/config.txt

module load Miniconda3
conda activate sqanti3


## Set input/output files
isoform_gff="${TRANSCRIPTOMEDIR}/Isoforms/collapsed.gff"
flnc_count="${TRANSCRIPTOMEDIR}/Isoforms/collapsed.flnc_count.txt"

sqanti_outdir="${TRANSCRIPTOMEDIR}/Annotation/SQANTI_classification/"
mkdir -p "${sqanti_outdir}"


## Check if SQANTI script exists 
echo "Checking if sqanti3_qc.py script exists"
if [ ! -s "${SQANTI3}" ]; then
    echo "ERROR: SQANTI3 script not found: ${SQANTI3}"
    exit 1
fi


## Check if SQANTI3 classification has already been completed
if [ -s "${sqanti_outdir}/isoforms_classification.txt" ]; then
    echo "SQANTI Isoform Classification output exists - skipping..."
else
    echo "Running Isoform Classification..."

    ## Check required input files exist
    if [ ! -s "${isoform_gff}" ]; then
        echo "ERROR: Missing or empty collapsed isoform GFF: ${isoform_gff}"
        exit 1
    fi

    if [ ! -s "${flnc_count}" ]; then
        echo "ERROR: Missing or empty FLNC count file: ${flnc_count}"
        exit 1
    fi

    if [ ! -s "${ANNOTATIONGTF}" ]; then
        echo "ERROR: Missing or empty reference annotation: ${ANNOTATIONGTF}"
        exit 1
    fi

    if [ ! -s "${REFGENOME}" ]; then
        echo "ERROR: Missing or empty reference genome: ${REFGENOME}"
        exit 1
    fi

    # Find STAR junction files
    shopt -s nullglob
    junction_files=( "${JUNCTIONDIR}"/*SJ.out.tab )
    shopt -u nullglob

    if [ ${#junction_files[@]} -eq 0 ]; then
        echo "ERROR: No STAR junction files (SJ.out.tab) found in ${JUNCTIONDIR}"
        exit 1
    fi

    ## Convert junction file paths into comma-separated list
    printf -v coverage_arg '%s,' "${junction_files[@]}"
    coverage_arg="${coverage_arg%,}"

    echo "Found ${#junction_files[@]} STAR junction file(s)."
    echo "Running SQANTI3 classification..."

    ## Run SQANTI3 classification
    python "${SQANTI3}" \
        --isoforms "${isoform_gff}" \
        --refGTF "${ANNOTATIONGTF}" \
        --refFasta "${REFGENOME}" \
        --fl_count "${flnc_count}" \
        --coverage "${coverage_arg}" \
        -d "${sqanti_outdir}" \
        --cpus "${SLURM_CPUS_PER_TASK}" \
        --include_ORF

    echo
    echo "Completed SQANTI3 classification..."

fi

echo "Master Transcriptome classification job completed"

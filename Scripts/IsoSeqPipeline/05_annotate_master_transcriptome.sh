#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=01:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/ClassifyIsoforms-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/ClassifyIsoforms-%j.err
#SBATCH --job-name=ClassifyIsoforms

## bash script to automate classification of collapsed isoforms 
## script also prepares (sort and index) files required for isoform classification
    ## these include reference annotations and isoforms GTF/GFF, reference sequence FASTA
## Utilised tools from PacBio Pigeon Transcript Tool

## do not store any sensitive data use config file to specify filepaths etc.

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/IsoSeqPipeline/05_annotate_master_transcriptome.sh

set -euo pipefail

echo "Starting Annotate Master Transcriptome job"
echo ""

## Load config and environment 
source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools 


## Set output dir
mkdir -p "${TRANSCRIPTOMEDIR}/Annotation/" 

out_annotate="${TRANSCRIPTOMEDIR}/Annotation/collapsed_classification.txt"
out_filtered="${out_annotate%.txt}.filtered_lite_classification.txt"
out_saturation="${TRANSCRIPTOMEDIR}/Annotation/saturation.txt"


## Set input file paths 
sorted_ref_gtf="${ANNOTATIONGTF%.gtf}.sorted.gtf"

isoform_gff="${TRANSCRIPTOMEDIR}/Isoforms/collapsed.gff"
sorted_isoform_gff="${isoform_gff%.gff}.sorted.gff"

flnc_count="${TRANSCRIPTOMEDIR}/Isoforms/collapsed.flnc_count.txt"


## Step 1: Prepare (sort and index) input files (if not done already)
# Step 1a: Prepare reference files 
if [ -s "${sorted_ref_gtf}" ]; then 
    echo "Reference GTF input already prepared - skipping..."
else
    echo "Preparing Reference GTF input file..."
    pigeon prepare ${ANNOTATIONGTF} ${REFGENOME} 
fi

# Validate output
if [ ! -s "${sorted_ref_gtf}" ]; then
    echo "ERROR: Reference GTF preparation failed"
    exit 1
fi


# Step 1b:  Prepare transcript isoforms GFF (output from isoform collapse)  
if [ -s "${sorted_isoform_gff}" ]; then 
    echo "Collapsed Isoform GFF input already prepared - skipping..."
else
    echo "Preparing Collapsed Isoform GFF input file..."
    pigeon prepare ${TRANSCRIPTOMEDIR}/Isoforms/collapsed.gff
fi

# Validate output 
if [ ! -s "${sorted_isoform_gff}" ]; then
    echo "ERROR: Isoform GFF preparation failed"
    exit 1
fi


## Step 2: Classify Isoforms into categories 
# if valid output file exists, skip step  
if [ -s "${out_annotate}" ]; then 
    echo "Isoform Classification output exists - skipping..."
else
    echo "Running Isoform Classification..."

    # Check if required input file exists
    if [ ! -s "${flnc_count}" ]; then
        echo "ERROR: Missing or empty FLNC count file in ${TRANSCRIPTOMEDIR}/Isoforms"
        exit 1
    fi

    # Run Pigeon Classify 
    pigeon classify \
        --fl ${flnc_count} \
        --out-dir ${TRANSCRIPTOMEDIR}/Annotation/ \
        --num-threads ${SLURM_CPUS_PER_TASK:-8} \
        ${sorted_isoform_gff} \
        ${sorted_ref_gtf} \
        ${REFGENOME} 
        
fi


## Step 3: Filter isoforms from the classification output
if [ -s "${out_filtered}" ]; then 
    echo "Filtered Isoform Classification file exists - skipping..."
else
    echo "Filtering Isoform classification file..."
    pigeon filter ${out_annotate} --isoforms ${sorted_isoform_gff}
fi


## Step 4: Report gene saturation [optional]
if [ -s "${out_saturation}" ]; then 
    echo "Gene saturation report exists - skipping..."
else
    echo "Generating gene saturation report..."
    pigeon report --exclude-singletons ${out_filtered} ${out_saturation}
fi


echo 
echo "Completed Transcriptome Annotation job"

## End of script 

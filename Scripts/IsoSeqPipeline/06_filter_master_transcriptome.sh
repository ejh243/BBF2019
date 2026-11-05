#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=02:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/FilterIsoforms-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/FilterIsoforms-%j.err
#SBATCH --job-name=FilterIsoforms

## bash script to automate filtering of classified isoforms 
## script runs pigeon pipeline and patches junctions statstics from SQANTI classification before filtering

## do not store any sensitive data use config file to specify filepaths etc.

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/IsoSeqPipeline/06_filter_master_transcriptome.sh

set -euo pipefail

echo "Starting Filter Master Transcriptome job"
echo ""

## Load config and environment 
source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools 

py_script="${SCRIPTSDIR}/IsoSeqPipeline/patch_coverage_data_into_pigeon.py"


## Set output directories
sqanti_outdir="${TRANSCRIPTOMEDIR}/Annotation/SQANTI_classification"
pigeon_outdir="${TRANSCRIPTOMEDIR}/Annotation/Pigeon_classification"
final_outdir="${TRANSCRIPTOMEDIR}/Annotation/Final_isoform_classification"

mkdir -p "${pigeon_outdir}" "${final_outdir}"

## Set output filenames 
sqanti_classify="${sqanti_outdir}/isoforms_classification.txt"
pigeon_classify="${pigeon_outdir}/collapsed_classification.txt"

out_patched="${final_outdir}/patched_classification.txt"
out_filtered="${out_patched%.txt}.filtered_lite_classification.txt"


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
    pigeon prepare "${TRANSCRIPTOMEDIR}/Isoforms/collapsed.gff"
fi
# Validate output 
if [ ! -s "${sorted_isoform_gff}" ]; then
    echo "ERROR: Isoform GFF preparation failed"
    exit 1
fi


## Step 2: Classify Isoform structures using Pigeon  
# if valid output file exists, skip step  
if [ -s "${pigeon_classify}" ]; then 
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
        --fl "${flnc_count}" \
        --out-dir "${pigeon_outdir}/" \
        --num-threads ${SLURM_CPUS_PER_TASK:-8} \
        "${sorted_isoform_gff}" \
        "${sorted_ref_gtf}" \
        "${REFGENOME}" 
fi
# Validate Pigeon classification output
if [ ! -s "${pigeon_classify}" ]; then
    echo "ERROR: Pigeon classification output missing or empty"
    exit 1
fi


## Step 3: Run py script to patch junction coverage data into Pigeon 
echo "Patching SQANTI3 junction coverage statistics into Pigeon classification..."
python "${py_script}" \
    --pigeon-classify "${pigeon_classify}" \
    --sqanti-classify "${sqanti_classify}" \
    --output "${out_patched}"

# Move Pigeon's junctions.txt into final_outdir to aid with filtering 
cp "${pigeon_outdir}/collapsed_junctions.txt" "${final_outdir}/patched_junctions.txt" 


## Step 4: Filter isoforms from the classification output
# Excluding any filters for monoexonic transcripts 
if [ -s "${out_filtered}" ]; then 
    echo "Filtered Isoform Classification file exists - skipping..."

else
    # Run filter 
    echo "Filtering Isoform classification file..."
    pigeon filter ${out_patched} --isoforms "${sorted_isoform_gff}"

    # Validate filtered output
    if [ ! -s "${out_filtered}" ]; then
        echo "ERROR: Pigeon filtering failed or output file missing"
        exit 1
    fi

    # Generate saturation report 
    echo "Generating gene/isoform saturation report..."
    pigeon report "${out_filtered}" "${final_outdir}/saturation_report.txt"
    pigeon report --exclude-singletons "${out_filtered}" "${final_outdir}/saturation_report_wo_singletons.txt"
fi


echo 
echo "Completed Transcriptome Annotation job"

## End of script 

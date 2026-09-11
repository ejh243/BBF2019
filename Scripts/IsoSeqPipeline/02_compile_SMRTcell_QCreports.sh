#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=00:10:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/SMRTcell_QC-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/SMRTcell_QC-%j.err
#SBATCH --job-name=SMRTcell_QC

## bash script to automate compilation of per SMRT cell QC metrics 
## runs the parse_SMRTcell_QC_metrics.py script 
## requires the following files to be present: 
    ## ${PROCESSEDDIR}/CCS/<subreads>.ccs_report.txt and 
    ## ${PROCESSEDDIR}/Refine/<subreads>.flnc.filter_summary.report.json 
    
## do not store any sensitive data use config file to specify filepaths etc.

## this script needs to be submitted from the main repository folder
## Usage: sbatch Scripts/IsoSeqPipeline/02_compile_SMRTcell_QCreports.sh

set -euo pipefail

echo "Starting Iso-Seq QC compilation" 
echo ""

## Load required software 
source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools  


## Set relevant filepaths 
py_script="${SCRIPTSDIR}/IsoSeqPipeline/parse_SMRTcell_QC_metrics.py"

ccs_dir="${PROCESSEDDIR}/CCS"
refine_dir="${PROCESSEDDIR}/Refine"

qc_dir="${PROCESSEDDIR}/QC"
mkdir -p ${qc_dir}

OUT_TSV="${qc_dir}/SMRTcell_QC_summary.tsv" 


# Write column headers on the output TSV (once)
echo -e \
"SMRT_cell_ID\t\
ZMWs_input\t\
ZMWs_passing_CCS\t\
ZMWs_passing_CCS(%)\t\
ZMWs_failing_CCS\t\
ZMWs_failing_CCS(%)\t\
Failed_CCS_lacking_full_passes\t\
Failed_CCS_lacking_full_passes(%)\t\
Failed_CCS_below_min_RQ\t\
Failed_CCS_below_min_RQ(%)\t\
Failed_coverage_drop\t\
Failed_coverage_drop(%)\t\
Failed_draft_generation\t\
Failed_draft_generation(%)\t\
FL_reads\t\
FLNC_reads\t\
FLNC_rate_(%)\t\
FLNC+polyA_reads\t\
FLNC+polyA_(%)" \
> "${OUT_TSV}"


# Counting number of SMRTcells found 
ccs_bams=("${ccs_dir}"/*.ccs.bam)

total_cells=${#ccs_bams[@]}
processed_cells=0
missing_ccs=0
missing_refine=0


# Loop over SMRT cells (derived from CCS.BAM file)
echo "Found ${total_cells} SMRT cells (from processed CCS BAM files)"
echo ""

for sample in "${ccs_bams[@]}"; do

    # SMRT cell ID is defined by CCS BAM
    smrt_cell_id=$(basename "${sample}" .ccs.bam)
    echo "Processing SMRT cell: ${smrt_cell_id}" 

    # Expected derived files
    ccs_report="${ccs_dir}/${smrt_cell_id}.ccs_report.txt"
    refine_report="${refine_dir}/${smrt_cell_id}.flnc.filter_summary.report.json"

    # Check if CCS report is present 
    if [[ ! -f "${ccs_report}" ]]; then
        echo "  Missing Isoseq CCS report for ${smrt_cell_id}, skipping." 
        missing_ccs=$((missing_ccs + 1))
        continue
    fi

    # Check if Refine report is present
    if [[ ! -f "${refine_report}" ]]; then
        echo "  Missing Isoseq Refine report for ${smrt_cell_id}, skipping." 
        missing_refine=$((missing_refine + 1))
        continue
    fi

    # Run python script to compile relevant SMRTcell QC metrics 
    python "${py_script}" \
        --smrt-cell "${smrt_cell_id}" \
        --ccs-report "${ccs_report}" \
        --refine-report "${refine_report}" \
        >> "${OUT_TSV}"

    processed_cells=$((processed_cells + 1))

done

echo ""
echo "QC compilation complete" 

echo ""
echo "Iso-Seq QC compilation summary"
echo "----------------------------------------"
echo "Total SMRT cells detected : ${total_cells}"
echo "Successfully processed    : ${processed_cells}"
echo "Missing CCS reports       : ${missing_ccs}"
echo "Missing refine reports    : ${missing_refine}"
echo "Output QC table           : ${OUT_TSV}"

# End of script 
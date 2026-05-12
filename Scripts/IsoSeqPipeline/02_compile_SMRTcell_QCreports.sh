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
#SBATCH --output=/lustre/home/vs455/LogFiles/SMRTcell_QC-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/SMRTcell_QC-%j.err
#SBATCH --job-name=SMRTcell_QC

## bash script to automate compilation of per SMRT cell QC metrics 
## requires the following files to be present: 
    ## ${PROCESSEDDIR}/CCS/<subreads>.ccs_report.txt and 
    ## ${PROCESSEDDIR}/Refine/<subreads>.flnc.filter_summary.report.json 
## do not store any sensitive data use config file to specify filepaths etc.

# Usage: sbatch Scripts/IsoSeqPipeline/02_compile_SMRTcell_QCreports.sh

source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools  


py_script="${SCRIPTSDIR}/IsoSeqPipeline/parse_SMRTcell_QC_metrics.py"

ccs_dir="${PROCESSEDDIR}/CCS"
refine_dir="${PROCESSEDDIR}/Refine"

OUT_TSV="${PROCESSEDDIR}/IsoSeq_QC_SMRTcells.tsv" 


# Write header (once)
echo -e \
"SMRT_cell_ID\tZMWs_input\tZMWs_passing_CCS\tZMWs_passing_CCS(%)\t\
ZMWs_failing_CCS\tZMWs_failing_CCS(%)\t\
Failed_CCS_lacking_full_passes\tFailed_CCS_lacking_full_passes(%)\t\
Failed_CCS_below_min_RQ\tFailed_CCS_below_min_RQ(%)\t\
FL_reads\tFLNC_reads\tFLNC_rate_(%)\t\
FLNC+polyA_reads\tFLNC+polyA_(%_of_input)" \
> "${OUT_TSV}"

# ----------------------------------------
# Loop over SMRT cells (derived from CCS reports)
# ----------------------------------------

for sample in "${ccs_dir}"/*.ccs.bam; do

    # SMRT cell ID is defined by CCS BAM
    smrt_cell_id=$(basename "${sample}" .ccs.bam)

    # Expected derived files
    ccs_report="${ccs_dir}/${smrt_cell_id}.ccs_report.txt"
    refine_report="${refine_dir}/${smrt_cell_id}.flnc.filter_summary.report.json"

    # Sanity checks
    if [[ ! -f "${ccs_report}" ]]; then
        echo "Missing Isoseq CCS report for ${smrt_cell_id}, skipping." >&2
        continue
    fi

    if [[ ! -f "${refine_report}" ]]; then
        echo "Missing Isoseq Refine report for ${smrt_cell_id}, skipping." >&2
        continue
    fi

    echo "Processing SMRT cell: ${smrt_cell_id}" >&2

    python "${py_script}" \
        --smrt-cell "${smrt_cell_id}" \
        --ccs-report "${ccs_report}" \
        --refine-report "${refine_report}" \
        >> "${OUT_TSV}"

done


echo "QC compilation complete: ${OUT_TSV}" >&2

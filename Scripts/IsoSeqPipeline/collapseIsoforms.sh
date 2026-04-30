#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=05:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/Collapse_test-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/Collapse_test-%j.err
#SBATCH --job-name=Collapse_test2

## bash script to automate collapsing of mapped transcripts into unique isoforms
## do not store any sensitive data use config file to specify filepaths etc.

# Note: The optional <flnc.bam> input is required to get the correct FLNC counts 
#       for bulk Iso-Seq in the flnc_count.txt supplemental file.

# Usage: sbatch Scripts/IsoSeqPipeline/collapseIsoforms.sh

source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools  


# Running with default isoform collapse logic (less strict)
if [ ! -f ${MASTERTRANSCRIPTOME}/collapsed_default.gff ]; then
    echo "Running isoform collapse with newer settings"
   
    isoseq collapse \
        --do-not-collapse-extra-5exons \
        --max-5p-diff 50 \
        --max-3p-diff 100 \
        --num-threads ${SLURM_CPUS_PER_TASK} \
        ${ALIGNEDDIR}/mapped.bam \
        ${PROCESSEDDIR}/Refine/flnc.fofn \
        ${MASTERTRANSCRIPTOME}/collapsed_default.gff

else
    echo "Collapse output file exists - skipping"
fi


# Running with legacy isoform collapse logic (stricter)  
# Can be removed later 
if [ ! -f ${MASTERTRANSCRIPTOME}/collapsed_legacy.gff ]; then
    echo "Running isoform collapse with legacy settings"
   
    isoseq collapse \
        --do-not-collapse-extra-5exons \
        --max-5p-diff 5 \
        --max-3p-diff 5 \
        --num-threads ${SLURM_CPUS_PER_TASK} \
        ${ALIGNEDDIR}/mapped.bam \
        ${PROCESSEDDIR}/Refine/flnc.fofn \
        ${MASTERTRANSCRIPTOME}/collapsed_legacy.gff

else
    echo "Legacy collapse output file exists - skipping"
fi





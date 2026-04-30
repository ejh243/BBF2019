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
#SBATCH --output=/lustre/home/vs455/LogFiles/Pigeon_classify_test-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/Pigeon_classify_test-%j.err
#SBATCH --job-name=Pigeon_classify_test

## bash script to classify isoform into categories 
## do not store any sensitive data use config file to specify filepaths etc.

# Usage: sbatch Scripts/PigeonPipeline/classifyIsoforms.sh

source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools  


# Run Isoform Classification step 
pigeon classify \
    ${MASTERTRANSCRIPTOME}/collapsed_default.sorted.gff \
    ${RESOURCESDIR}/gencode.v38.annotation.sorted.gtf \
    ${REFGENOME} \
    --fl ${MASTERTRANSCRIPTOME}/collapsed_default.flnc_count.txt
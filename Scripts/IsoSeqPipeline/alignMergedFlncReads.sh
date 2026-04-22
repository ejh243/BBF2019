#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/Alignment_test-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/Alignment_test-%j.err
#SBATCH --job-name=Alignment_test

## bash script to automate alignment of isoseq data
## assumes all SMRT cell runs have been merged before clustering
## do not store any sensitive data use config file to specify filepaths etc.

# Usage: sbatch Scripts/IsoSeqPipeline/alignMergedFlncReads.sh

source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools  

mkdir -p ${ALIGNEDDIR} # create dir if not present 


# Run PacBio's minimap2 alignment to map all FLNC reads to reference 

if [ ! -f ${ALIGNEDDIR}/flnc.fofn.mapped.bam ]; then
    echo "Running alignment of all FLNC reads"

    pbmm2 align \
        --preset ISOSEQ \
        --sort \
        --num-threads ${SLURM_CPUS_PER_TASK} \
        --log-file ${ALIGNEDDIR}/alignment.log \
        ${REFGENOME} \
        ${PROCESSEDDIR}/Cluster2/clustered_flnc_fofn.bam \
        ${ALIGNEDDIR}/flnc_fofn_mapped.bam 

else
    echo "Alignment file exists - skipping"
fi

echo "Done"

# End of script
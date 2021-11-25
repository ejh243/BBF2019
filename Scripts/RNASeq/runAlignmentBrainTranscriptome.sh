#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/alignSRBrain-%A_%a.o
#SBATCH --error=LogFiles/alignSRBrain-%A_%a.e
#SBATCH --job-name=alignSRBrain-%A_%a.e
#SBATCH --array=0-49%10 ## runs 50 jobs with 10 at any one time


module load STAR

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file
echo "Loading config file: "
source ./Config/config.txt

cd ${RNASeqDIR}

FQFILES=($(find . -name '*[rR]1*q.gz' -not -path "./Trimmed/*"))

echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	

sample=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleName=$(basename ${sample%_[GCTA]*})

cd ${SCRIPTSDIR}

## align to MergedTranscriptome
#sh RNASeq/alignShortReadCortexTranscriptome.sh ${sampleName}

module load RSEM
## gene & isofrom counts for to MergedTranscriptome
sh RNASeq/rsemBrainTranscriptome.sh ${sampleName}

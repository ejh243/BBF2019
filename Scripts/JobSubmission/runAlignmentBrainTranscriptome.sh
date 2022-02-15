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

## load config file
echo "Loading config file: "
source ./Config/config.txt

RNASEQDIR=$1
PROJECT=$2
GENECOUNTDIR=${GENECOUNTPATH}/${PROJECT}

mkdir -p ${GENECOUNTDIR}

FQFILES=($(find ${RNASEQDIR} -maxdepth 1 -name '*[rR]1*q.gz' ))
sample=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleName=$(basename ${sample%[rR]1*})

module load STAR
module load RSEM
## gene and isoform counts for GENCODE transcripts
sh Scripts/RNASeq/rsemBrainTranscriptome.sh ${sampleName} ${RNASEQDIR} ${GENECOUNTDIR}

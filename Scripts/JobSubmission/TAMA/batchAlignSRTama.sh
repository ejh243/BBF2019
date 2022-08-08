#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/alignSRTama.o
#SBATCH --error=LogFiles/alignSRTama.e
#SBATCH --job-name=alignSRTama
#SBATCH --array=0-34%10 ## runs multiple jobs with 10 at any one time

# this script needs to be submitted from the main repository folder



## load config file
echo "Loading config file: "
source ./Config/config.txt

RNASEQDIR=$1
PROJECT=$2
ALIGNEDDIR=${ALIGNEDPATH}/${PROJECT}
GENECOUNTDIR=${GENECOUNTPATH}/${PROJECT}

mkdir -p ${GENECOUNTDIR}

FQFILES=($(find ${RNASEQDIR} -maxdepth 1 -name '*[rR]1*q.gz' ))
sample=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleName=$(basename ${sample%[rR]1*})

module load STAR
module load RSEM
## align with brain transcriptome
sh Scripts/RNASeq/alignShortReadBrainTranscriptome.sh ${sampleName} ${RNASEQDIR} ${ALIGNEDDIR} ${GENECOUNTDIR} BrainTAMAPre


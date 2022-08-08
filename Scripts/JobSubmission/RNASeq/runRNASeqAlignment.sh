#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/alignShortRead-%A_%a.o
#SBATCH --error=LogFiles/alignShortRead-%A_%a.e
#SBATCH --job-name=alignShortRead-%A_%a.e


module load STAR

## load config file provided on command line when submitting job
echo "Loading config file: "
source ./Config/config.txt

RNASEQDIR=$1
PROJECT=$2
ALIGNEDDIR=${ALIGNEDPATH}/${PROJECT}
GENECOUNTDIR=${GENECOUNTPATH}/${PROJECT}
QCDIR=${ALIGNEDDIR}/rnaseqc

mkdir -p ${ALIGNEDDIR}
mkdir -p ${GENECOUNTDIR}
mkdir -p ${QCDIR}

FQFILES=($(find ${RNASEQDIR} -maxdepth 1 -name '*[rR]1*q.gz' ))

echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	

sample=${FQFILES[${SLURM_ARRAY_TASK_ID}]}
sampleName=$(basename ${sample%[rR]1*})

## align to gencode
sh Scripts/RNASeq/alignShortReadData.sh ${sampleName} ${RNASEQDIR} ${ALIGNEDDIR}

module load RSEM
## gene and isoform counts for GENCODE transcripts
sh Scripts/RNASeq/rsemGENCODE.sh ${sampleName} ${RNASEQDIR} ${GENECOUNTDIR}

## 
module load Miniconda2
source ./Config/config.txt
source activate rnaseqc

sh Scripts/RNASeq/rnaseqQC.sh ${sampleName} ${ALIGNEDDIR} ${QCDIR}
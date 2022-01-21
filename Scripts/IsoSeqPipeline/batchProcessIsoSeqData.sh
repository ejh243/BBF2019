#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/PreprocessIsoseq3s-%A_%a.o
#SBATCH --error=LogFiles/PreprocessIsoseq3s-%A_%a.e
#SBATCH --job-name=PreprocessIsoseq3s-%A_%a.e
#SBATCH --array=0-32%10 ## runs multiple jobs with 10 at any one time

# this script needs to be submitted from the main repository folder

source ./Config/config.txt

echo "Changing Folder to: "
echo $DATADIR
cd $DATADIR/

## this command can be used to process all relevant files in DATADIR 
## if only a subset need to be processed provide list in FilesToProcess.txt and hash out line below.

samples=($(ls *.subreads.bam))

echo "Samples to process: " ${#samples[@]}

sample=${samples[${SLURM_ARRAY_TASK_ID}]}


## run first steps on each smrt cell individually
mkdir -p ${PROCESSEDDIR}/CCS
mkdir -p ${PROCESSEDDIR}/Lima
mkdir -p ${PROCESSEDDIR}/Refine
mkdir -p ${PROCESSEDDIR}/Cluster
mkdir -p ${PROCESSEDDIR}/Polish


echo "Changing Folder to: "
echo ${SCRIPTSDIR}/IsoSeqPipeline
cd ${SCRIPTSDIR}/IsoSeqPipeline


module load Miniconda2
source activate isoseq

## output version of ccs
ccs --version
## output version of lima
lima --version
#sh ./processIsoSeqSMRTcells.sh ${sample}

module load minimap2
#sh ./alignIsoSeqSMRTcells.sh ${sample}


module purge
module load Miniconda2
source activate anaCogent

mkdir -p ${ALIGNEDDIR}/Collapsed/


#sh ./filterIsoSeqSMRTcells.sh ${sample}
module purge
module load STAR
module load RSEM

mkdir -p ${RSEMREFDIR}
mkdir -p ${GENECOUNTSDIR}/RSEM/PersonalTranscriptome/

sh ./alignShortReadSMRTcells.sh ${sample}

module purge
module load Miniconda2
source activate SQANTI3.env

sh ./sqanti3QCIsoSeqSMRTcells.sh ${sample}

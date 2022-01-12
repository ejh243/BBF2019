#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p sq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/PreprocessIsoseq3s-%A_%a.o
#SBATCH --error=LogFiles/PreprocessIsoseq3s-%A_%a.e
#SBATCH --job-name=PreprocessIsoseq3s-%A_%a.e
#SBATCH --array=19-32%10 ## runs multiple jobs with 10 at any one time

echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

module load Miniconda2
source activate isoseq
source ./Config/config.txt

## output version of ccs
ccs --version
## output version of lima
lima --version

echo "Changing Folder to: "
echo $DATADIR


cd $DATADIR/


## this command can be used to process all relevant files in DATADIR 
## if only a subset need to be processed provide list in FilesToProcess.txt and hash out line below.

samples=($(ls *.subreads.bam))

echo "Samples to process: " ${#samples[@]}

sample=${samples[${SLURM_ARRAY_TASK_ID}]}


## run first steps on each sample individually
mkdir -p ${PROCESSEDDIR}/CCS
mkdir -p ${PROCESSEDDIR}/Lima
mkdir -p ${PROCESSEDDIR}/Refine

echo "Changing Folder to: "
echo ${SCRIPTSDIR}/IsoSeqPipeline
cd ${SCRIPTSDIR}/IsoSeqPipeline
sh ./ProcessIsoSeqFiles.sh ${sample}
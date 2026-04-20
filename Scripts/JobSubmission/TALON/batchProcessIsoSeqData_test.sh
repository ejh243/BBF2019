#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # enter email address
#SBATCH --output=/lustre/home/vs455/LogFiles/PreprocessIsoseq_test-%A_%a.out 
#SBATCH --error=/lustre/home/vs455/LogFiles/PreprocessIsoseq_test-%A_%a.err 
#SBATCH --job-name=PreprocessIsoseq_test
#SBATCH --time=05:00:00
#SBATCH --array=0-1%2 ## runs multiple jobs with 10 at any one time (set 0-34%10 for full cohort)

# this script needs to be submitted from the main repository folder
# Example Usage: sbatch Scripts/JobSubmission/TALON/batchProcessIsoSeqData_test.sh


source ./Config/config.txt

echo "Changing Folder to: "
echo $DATADIR
cd $DATADIR/

## this command can be used to process all relevant files in DATADIR 
## if only a subset need to be processed provide list in FilesToProcess.txt and hash out line below.

samples=( *.subreads.bam )

echo "Samples to process: " ${#samples[@]}

sample=${DATADIR}/${samples[${SLURM_ARRAY_TASK_ID}]}


## run initial steps on each SMRT cell individually
mkdir -p ${PROCESSEDDIR}/CCS
mkdir -p ${PROCESSEDDIR}/Lima
mkdir -p ${PROCESSEDDIR}/Refine


echo "Changing Folder to: "
echo ${SCRIPTSDIR}/IsoSeqPipeline
cd ${SCRIPTSDIR}/IsoSeqPipeline


module load Miniconda3
source activate isoseq_tools   

echo "software tools used"
## output version of isoseq
isoseq --version
## output version of ccs
ccs --version
## output version of lima
lima --version

# Run the Iso Seq per-SMRT cell pre-processing script 
bash ./processIsoSeqSMRTcells_updated.sh ${sample} 


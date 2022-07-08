#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/tamaProcess-%A_%a.o
#SBATCH --error=LogFiles/tamaProcess-%A_%a.e
#SBATCH --job-name=tamaProcess-%A_%a.e
#SBATCH --array=0-34%10 ## runs multiple jobs with 10 at any one time

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

echo "Changing Folder to: "
echo ${SCRIPTSDIR}/tamaPipeline
cd ${SCRIPTSDIR}/tamaPipeline

module load minimap2
module load BEDTools
module load SAMtools

mkdir -p ${ALIGNEDDIR}/TAMA/

sh ./alignFLNC.sh ${sample}

module purge
module load Biopython/1.72-foss-2018b-Python-2.7.15
module load Pysam/0.15.1-foss-2018b-Python-2.7.15
module load SAMtools

mkdir -p ${ALIGNEDDIR}/TAMA/Collapsed/

sh ./tamaCollapse.sh ${sample}
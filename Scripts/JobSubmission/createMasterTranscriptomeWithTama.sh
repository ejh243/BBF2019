#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/tamaMerge.o
#SBATCH --error=LogFiles/tamaMerge.e
#SBATCH --job-name=tamaMerge

# this script needs to be submitted from the main repository folder

source ./Config/config.txt
module purge
module load Biopython/1.72-foss-2018b-Python-2.7.15
module load Pysam/0.15.1-foss-2018b-Python-2.7.15
module load SAMtools

mkdir -p ${MASTERTRANSCRIPTOME}/TAMA

cd ${ALIGNEDDIR}/TAMA/Collapsed/
python ${SOFTWAREPATH}/tama/tama_merge.py -f smrtcells.txt -p ${MASTERTRANSCRIPTOME}/TAMA/pfc_merge_smrt_all -d merge_dup -a 100 -z 100
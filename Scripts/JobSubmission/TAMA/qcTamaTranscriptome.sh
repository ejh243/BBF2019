#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/QCtamaTranscriptome
#SBATCH --error=LogFiles/QCtamaTranscriptome
#SBATCH --job-name=QCtamaTranscriptome


module purge
module load Miniconda2
source activate SQANTI3.env

cd ${SCRIPTSDIR}/tamaPipeline

sh sqantiqc.sh



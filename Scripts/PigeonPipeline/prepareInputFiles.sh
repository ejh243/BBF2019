#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=05:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks=1 # specify number of tasks per node
#SBATCH --cpus-per-task=16 # full node utilisation 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=v.suresh@exeter.ac.uk # email me at job completion
#SBATCH --output=/lustre/home/vs455/LogFiles/Pigeon_prepare_test-%j.out
#SBATCH --error=/lustre/home/vs455/LogFiles/Pigeon_prepare_test-%j.err
#SBATCH --job-name=Pigeon_prepare_test

## bash script to prepare (sort and index) files required for isoform classification
## these include reference annotations and isoforms GTF/GFF, reference sequence FASTA
## multiple files may be listed on the command line, or a file of filenames may be provided.
## do not store any sensitive data use config file to specify filepaths etc.

# Usage: sbatch Scripts/PigeonPipeline/prepareInputFiles.sh

source ./Config/config.txt

module load Miniconda3
source activate isoseq_tools  

# Create and move to Annotation directory  
mkdir -p ${ANNOTATIONDIR} 
cd ${ANNOTATIONDIR} 


# Prepare reference files 
pigeon prepare --log-file ${ANNOTATIONDIR}/prepare_ref.log ${GENCODEGTF} ${REFGENOME} 


# Prepare transcript isoforms GFF (output from isoform collapse)  
pigeon prepare --log-file ${ANNOTATIONDIR}/prepare_iso_gff.log ${MASTERTRANSCRIPTOME}/collapsed.gff




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
module load BEDTools

mkdir -p ${MASTERTRANSCRIPTOME}/TAMA

cd ${ALIGNEDDIR}/TAMA/Collapsed/
python ${SOFTWAREPATH}/tama/tama_merge.py -f smrtcells.txt -p ${MASTERTRANSCRIPTOME}/TAMA/pfc_merge_smrt_all -d merge_dup -a 100 -z 100

# calculate read 
python ${SOFTWAREPATH}/tama/tama_go/read_support/tama_read_support_levels.py -f ${MASTERTRANSCRIPTOME}/TAMA/readSupport.txt -o ${MASTERTRANSCRIPTOME}/TAMA/pfc_merge_smrt_all_counts -m ${MASTERTRANSCRIPTOME}/TAMA/pfc_merge_smrt_all_merge.txt


# filter singletons
cd ${MASTERTRANSCRIPTOME}/TAMA/
python ${SOFTWAREPATH}/tama/tama_go/filter_transcript_models/tama_remove_single_read_models_levels.py -b pfc_merge_smrt_all.bed -r pfc_merge_smrt_all_counts_read_support.txt -o pfc_merge_smrt_all_nRead2

# remove transcripts that may be fragments of longer transcripts
python ${SOFTWAREPATH}/tama/tama_go/filter_transcript_models/tama_remove_fragment_models.py -f pfc_merge_smrt_all_nRead2.bed -o pfc_merge_smrt_all_nRead2_filtFrag

# convert to gtf
python ${SOFTWAREPATH}/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py pfc_merge_smrt_all_nRead2_filtFrag.bed pfc_merge_smrt_all_nRead2_filtFrag.gtf


module purge
module load STAR
module load RSEM

sh ${SCRIPTSDIR}/tamaPipeline/alignShortRead.sh
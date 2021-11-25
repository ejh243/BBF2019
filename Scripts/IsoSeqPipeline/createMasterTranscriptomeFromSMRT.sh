#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the high memory queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/MasterTranscriptome.o
#SBATCH --error=LogFiles/MasterTranscriptome.e
#SBATCH --job-name=MasterTranscriptome

cd $SLURM_SUBMIT_DIR

module load minimap2
module load SAMtools

source ./Config/config.txt

## this scripts combines qc'd isoform data across SMRT  cells using TALON
## recommended for moderate-to-large (ex: >4) sample merging where a reference transcriptome readily exists

	 
# making MD-tagged SAM files for each sample
FASTAFILES=($(ls  ${ALIGNEDDIR}/Collapsed/m*/SQANTI3/*_classification.filtered_lite.fasta))

mkdir -p ${MASTERTRANSCRIPTOME}/TALON
cd ${MASTERTRANSCRIPTOME}/TALON
mkdir -p InputSAM

## create config files
echo -n InputSAM/config.smrt.csv

for fasta in ${FASTAFILES[@]}
do
	sampleName=$(basename ${fasta})
	sampleName=${sampleName%_classification.filtered_lite.fasta}
	
	# align to ref genome
	minimap2 -ax splice -uf --secondary=no -C5 ${REFGENOME} -t30 \
		${fasta} > InputSAM/${sampleName}.collapsed.rep.fa.hg38.sam

	samtools calmd InputSAM/${sampleName}.collapsed.rep.fa.hg38.sam \
		${REFGENOME} \
		--output-fmt sam > InputSAM/${sampleName}.collapsed.rep.fa.hg38.MDtagged.sam
	
	echo ${sampleName},${sampleName},PacBio,InputSAM/${sampleName}.collapsed.rep.fa.hg38.MDtagged.sam >> InputSAM/config.smrt.csv
done





module purge
module load Miniconda2
source activate talon

# initializing the database

talon_initialize_database --f ${GENCODEGTF}\
     --g hg38 --a gencode38 --o pfc_merge_smrt

## add samples to the database
talon --f InputSAM/config.smrt.csv \
      --db pfc_merge_smrt.db \
      --build hg38 \
      --t 30 \
      --cov 0.95 \
      --identity 0.95 \
      --o pfc_merge_smrt

## summarise transcript numbers
talon_summarize \
       --db pfc_merge_smrt.db \
       --v \
       --o pfc_merge_smrt

## create an abundance matrix
talon_abundance \
       --db pfc_merge_smrt.db \
       -a gencode38 \
       --build hg38 \
       --o pfc_merge_smrt
	   
## create gtf
talon_create_GTF --db=pfc_merge_smrt.db \
                 --annot=gencode38 \
                 -b hg38 \
                 --observed \
                 --o pfc_merge_smrt



## compare transcriptome to that from Leung at al.
source Config/config.txt

BRAINGTF=${MASTERTRANSCRIPTOME}/TALON/pfc_merge_filter_talon_observedOnly.gtf

module purge
module load Miniconda2
source activate talon


gffcompare -o ${MASTERTRANSCRIPTOME}/TALON/ClassifyRefLeung -r ${REFGTF} ${BRAINGTF}



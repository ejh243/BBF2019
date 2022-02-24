
source $1 

## create STAR index for reference genome
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${STARIndex} \
--genomeFastaFiles ${REFGENOME} \
--sjdbGTFfile ${GENCODEGTF} \
--sjdbOverhang 99

## create STAR index for merged master transcriptome

BRAINGTF=${MASTERTRANSCRIPTOME}/TALON/pfc_merge_filter_talon_observedOnly.gtf
mkdir -p ${STARINDEXDIR}/Brain
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${STARINDEXDIR}/Brain \
--genomeFastaFiles ${REFGENOME} \
--sjdbGTFfile ${BRAINGTF} \
--sjdbOverhang 99


source $1 

## create STAR index for reference genome
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${STARIndex} \
--genomeFastaFiles ${REFGENOME} \
--sjdbGTFfile ${GENCODEGTF} \
--sjdbOverhang 99

## create STAR index for merged master transcriptome with TALON
BRAINGTF=${MASTERTRANSCRIPTOME}/TALON/pfc_merge_filter_talon_observedOnly.gtf
mkdir -p ${STARINDEXDIR}/BrainTALON
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${STARINDEXDIR}/BrainTALON \
--genomeFastaFiles ${REFGENOME} \
--sjdbGTFfile ${BRAINGTF} \
--sjdbOverhang 99   --limitSjdbInsertNsj 1091177

## create STAR index for QC'd merged master transcriptome with TAMA post SQANTIQC
BRAINGTF=${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/TAMAmerge_filter.filtered.gtf
mkdir -p ${STARINDEXDIR}/BrainTAMAPost
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${STARINDEXDIR}/BrainTAMAPost \
--genomeFastaFiles ${REFGENOME} \
--sjdbGTFfile ${BRAINGTF} \
--sjdbOverhang 99  --limitSjdbInsertNsj 1091177

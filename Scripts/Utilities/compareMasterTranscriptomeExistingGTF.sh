## compare transcriptome to that from Leung at al.
source Config/config.txt

BRAINGTF=${MASTERTRANSCRIPTOME}/TALON/pfc_merge_filter_talon_observedOnly.gtf

module purge
module load Miniconda2
source activate talon


gffcompare -o ${MASTERTRANSCRIPTOME}/TALON/CompareGTEX -r ${RESOURCESDIR}/GTEX/flair_filter_transcripts.gtf ${BRAINGTF}

gffcompare -o ${MASTERTRANSCRIPTOME}/TALON/CompareCHESS -r ${RESOURCESDIR}/CHESS/chess2.2.gff ${BRAINGTF}

gffcompare -o ${MASTERTRANSCRIPTOME}/TALON/CompareGENCODE -r ${GENCODEGTF} ${BRAINGTF}

BRAINGTF=${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/TAMAmerge_filter.filtered.gtf

gffcompare -o ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/CompareGTEX -r ${RESOURCESDIR}/GTEX/flair_filter_transcripts.gtf ${BRAINGTF}

gffcompare -o ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/CompareCHESS -r ${RESOURCESDIR}/CHESS/chess2.2.gff ${BRAINGTF}

gffcompare -o ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/CompareGENCODE -r ${GENCODEGTF} ${BRAINGTF}

gffcompare -o ${MASTERTRANSCRIPTOME}/CompareTAMAvsTALONPipeline -r ${BRAINGTF} ${MASTERTRANSCRIPTOME}/TALON/pfc_merge_filter_talon_observedOnly.gtf
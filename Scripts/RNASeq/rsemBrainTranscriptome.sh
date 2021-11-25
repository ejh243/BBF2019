sampleName=$1
BRAINGTF=${MASTERTRANSCRIPTOME}/TALON/pfc_merge_smrt_talon_observedOnly.gtf
	


cd ${RNASeqDIR}
TRIMFOLDER=Trimmed

star_f1=$(ls ${TRIMFOLDER}/${sampleName}_*[rR]1*f*z)
star_f2=$(ls ${TRIMFOLDER}/${sampleName}_*[rR]2*f*z)

nCPUS=$(($SLURM_CPUS_ON_NODE * $SLURM_CPUS_ON_NODE))

mkdir -p ${RSEMREFDIR}/talon_merged/talon_merged
#rsem-prepare-reference --gtf ${BRAINGTF} --star ${REFGENOME} ${RSEMREFDIR}/talon_merged/talon_merged

mkdir -p ${GENECOUNTSDIR}/RSEM/MergedTranscriptome/

rsem-calculate-expression --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/talon_merged/talon_merged ${GENECOUNTSDIR}/RSEM/MergedTranscriptome/talon_merged_lbb_${sampleName}


sampleName=$1


cd ${RNASeqDIR}
TRIMFOLDER=Trimmed

star_f1=$(ls ${TRIMFOLDER}/${sampleName}_*[rR]1*f*z)
star_f2=$(ls ${TRIMFOLDER}/${sampleName}_*[rR]2*f*z)

nCPUS=$(($SLURM_CPUS_ON_NODE * $SLURM_CPUS_ON_NODE))

#mkdir ${RSEMREFDIR}/GENCODEv38/
#rsem-prepare-reference --gtf ${GENCODEGTF} --star ${REFGENOME} ${RSEMREFDIR}/GENCODEv38/GENCODEv38

mkdir -p ${GENECOUNTSDIR}/RSEM/GENCODE/

rsem-calculate-expression --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/GENCODEv38/GENCODEv38 ${GENECOUNTSDIR}/RSEM/GENCODE/GENCODEv38_${sampleName}




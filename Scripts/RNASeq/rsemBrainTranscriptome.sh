
BRAINGTF=${MASTERTRANSCRIPTOME}/TALON/pfc_merge_smrt_talon_observedOnly.gtf

sampleName=$1
RNASEQDIR=$2
GENECOUNTSDIR=$3

TRIMDIR=trimmed

cd ${RNASEQDIR}

star_f1=$(ls ${TRIMDIR}/${sampleName}*[rR]1*f*z)
star_f2=$(ls ${TRIMDIR}/${sampleName}*[rR]2*f*z)

mkdir -p ${GENECOUNTSDIR}/RSEM/MergedTranscriptome/

nCPUS=$(($SLURM_CPUS_ON_NODE * $SLURM_CPUS_ON_NODE))

rsem-calculate-expression -p ${nCPUS} --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/talon_merged/talon_merged ${GENECOUNTSDIR}/RSEM/MergedTranscriptome/brainMasterGTF_${sampleName}

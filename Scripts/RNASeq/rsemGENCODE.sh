
sampleName=$1
RNASEQDIR=$2
GENECOUNTSDIR=$3

TRIMDIR=trimmed


cd ${RNASEQDIR}


star_f1=$(ls ${TRIMDIR}/${sampleName}*[rR]1*f*z)
star_f2=$(ls ${TRIMDIR}/${sampleName}*[rR]2*f*z)

nCPUS=$(($SLURM_CPUS_ON_NODE * $SLURM_CPUS_ON_NODE))


mkdir -p ${GENECOUNTSDIR}/RSEM/GENCODE/

rsem-calculate-expression --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/GENCODEv38/GENCODEv38 ${GENECOUNTSDIR}/RSEM/GENCODE/GENCODEv38_${sampleName}





sampleName=$1
RNASEQDIR=$2
GENECOUNTDIR=$3

TRIMDIR=trimmed


cd ${RNASEQDIR}


star_f1=$(ls ${TRIMDIR}/${sampleName}*[rR]1*f*z)
star_f2=$(ls ${TRIMDIR}/${sampleName}*[rR]2*f*z)

nCPUS=$(($SLURM_CPUS_ON_NODE * $SLURM_CPUS_ON_NODE))


mkdir -p ${GENECOUNTDIR}/RSEM/GENCODE/

if [ ! -f ${GENECOUNTDIR}/RSEM/GENCODE/GENCODEv38_${sampleName}.genes.results ]
then	
	echo "gene counting" ${sampleName}

	rsem-calculate-expression --star --star-gzipped-read-file --paired-end ${star_f1} ${star_f2} ${RSEMREFDIR}/GENCODEv38/GENCODEv38 ${GENECOUNTDIR}/RSEM/GENCODE/GENCODEv38_${sampleName}
fi



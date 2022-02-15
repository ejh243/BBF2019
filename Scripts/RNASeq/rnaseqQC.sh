
sampleName=$1
ALIGNEDRNA=$2
QCDIR=$3


cd ${SOFTWAREPATH}/rnaseqc

## run rnaseqc for each bam file
bamfile=${ALIGNEDRNA}/GENCODE_v38_${sampleName}Aligned.sortedByCoord.out.bam

${SOFTWAREPATH}/rnaseqc.v2.4.2.linux ${GENCODEGTF/%annotation.gtf/genes.gtf} ${bamfile} ${QCDIR} -s ${sampleName}

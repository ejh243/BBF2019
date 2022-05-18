
sampleName=$1
RNASEQDIR=$2
ALIGNEDRNA=$3

TRIMDIR=trimmed

##find trimmed files

cd ${RNASEQDIR}
star_f1=$(ls ${TRIMDIR}/${sampleName}*[rR]1*f*z)
star_f2=$(ls ${TRIMDIR}/${sampleName}*[rR]2*f*z)

mkdir -p ${ALIGNEDRNA}

module purge
module load STAR

if [ ! -f ${ALIGNEDRNA}/Brain_${sampleName}Log.final.out ]
then	
	echo "aligning" ${star_f1}

	## align with STAR using GENCODE v38 star index
	STAR --genomeDir ${STARINDEXDIR}/Brain --runThreadN 18 --readFilesIn ${star_f1},${star_f2} \
		--readFilesCommand zcat \
		--outFileNamePrefix ${ALIGNEDRNA}/Brain_${sampleName} \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes Standard
fi

	
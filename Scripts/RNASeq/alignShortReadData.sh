## process short read data and align to GENCODE

## generate QC report with fastqc
module load FastQC

sampleName=$1

cd ${RNASeqDIR}

f1=$(ls Sample_${sampleName}/raw_illumina_reads/${sampleName}_*[rR]1*f*z)
f2=$(ls Sample_*${sampleName}/raw_illumina_reads/${sampleName}_*[rR]2*f*z)

echo "found fastq files for sample" ${sampleName}
echo ${f1}
echo ${f2}


## trimmomatic
module purge
module load Trim_Galore/0.4.5-foss-2016b


TRIMFOLDER=Trimmed
mkdir -p ${TRIMFOLDER}


if [[ ! -f ${TRIMFOLDER}/${sampleName}_*[rR]1*_val_1*f*z ]] && [[ ! -f ${TRIMFOLDER}/${sampleName}_*[rR]2*_val_2*f*z ]]
then	
	echo "trimming" ${f1}
	## trim reads
	trim_galore --fastqc --gzip --suppress_warn -o ${TRIMFOLDER} --paired ${f1} ${f2}
fi

star_f1=$(ls ${TRIMFOLDER}/${sampleName}_*[rR]1*f*z)
star_f2=$(ls ${TRIMFOLDER}/${sampleName}_*[rR]2*f*z)

mkdir -p ${ALIGNEDRNA}

module purge
module load STAR

if [ ! -f ${ALIGNEDRNA}/GENCODE_v38_${sampleName}Log.final.out ]
then	
	echo "aligning" ${star_f1}

	## align with STAR using GENCODE v38 star index
	STAR --genomeDir ${STARIndex} --runThreadN 18 --readFilesIn ${star_f1},${star_f2} \
		--readFilesCommand zcat \
		--outFileNamePrefix ${ALIGNEDRNA}/GENCODE_v38_${sampleName} \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes Standard
fi

	
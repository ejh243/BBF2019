#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the high memory queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/ClusterIndividuals.o
#SBATCH --error=LogFiles/ClusterIndividuals.e
#SBATCH --job-name=ClusterIndividuals



echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

module load Miniconda2
source activate isoseq
source ./Config/config.txt

## create matched arrays of isoseq sample IDs, and brain ids where brain id is a merge of brainbank id and rna id
isoseqIDs=()
brainIDs=()
{
read ## to ignore first line
while IFS=',' read -ra ADDR
do
  isoseqIDs+=(${ADDR[2]%.subreads.bam})
  brainIDs+=(${ADDR[8]}"_"${ADDR[1]})
done } < "$SAMPLESHEET"

## create list of unique brainIDs
uniqIDs=($(printf "%s\n" "${brainIDs[@]}" | sort -u))

echo "Changing Folder to: "
echo $DATADIR

## SMRT cells will be merged for subsequent steps
cd $DATADIR/


## cluster data
mkdir -p ${PROCESSEDDIR}/Cluster
mkdir -p ${PROCESSEDDIR}/Polish



for sample in ${uniqIDs[@]}
do
	XMLFILES=()
	## need to create a fofn for this individual
	echo -n > ${PROCESSEDDIR}/Refine/${sample}_flnc.fofn
	for keys in "${!brainIDs[@]}"
	do 
		if [ ${brainIDs[$keys]} = ${sample} ];
		then 
			echo "Matched to "  ${isoseqIDs[$keys]}
			echo ${PROCESSEDDIR}/Refine/${isoseqIDs[$keys]}.flnc.bam >> ${PROCESSEDDIR}/Refine/${sample}_flnc.fofn
			XMLFILES+=(${isoseqIDs[$keys]}.subreadset.xml)
		fi
	done
	
	echo "Clustering " ${sample}	
	
	if [ ! -f ${PROCESSEDDIR}/Cluster/clustered_${sample}.bam ]
	then 
		isoseq3 cluster ${PROCESSEDDIR}/Refine/${sample}_flnc.fofn ${PROCESSEDDIR}/Cluster/clustered_${sample}.bam --verbose --use-qvs --log-file ${PROCESSEDDIR}/Cluster/merge_${sample}.log
	fi 
	
	echo "Polishing " ${sample}
	
	if [ ! -f ${PROCESSEDDIR}/Cluster/Polish/polished_${sample}.bam ]
	then

		# create merged xml
		dataset create --force --type SubreadSet --novalidate ${DATADIR}/${sample}.subreadset.xml ${XMLFILES[@]}

		isoseq3 polish ${PROCESSEDDIR}/Cluster/clustered_${sample}.bam ${sample}.subreadset.xml ${PROCESSEDDIR}/Polish/polished_${sample}.bam
	fi 

done



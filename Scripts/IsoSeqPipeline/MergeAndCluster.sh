#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the high memory queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/MergeAndCluster.o
#SBATCH --error=LogFiles/MergeAndCluster.e
#SBATCH --job-name=MergeAndCluster


echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

module load Miniconda2
source activate isoseq
source ./Config/config.txt

echo "Changing Folder to: "
echo $DATADIR

## SMRT cells will be merged for subsequent steps
cd $DATADIR/


## cluster data
mkdir -p ${PROCESSEDDIR}/Cluster
mkdir -p ${PROCESSEDDIR}/Polish

## first cluster each smrt cell
flnc=($(ls ${PROCESSEDDIR}/Refine/*.flnc.bam))
for file in ${flnc[@]}
do
	sample=$(basename ${file})
	sampleName=${sample%.flnc.bam}

	echo "Clustering " ${sampleName}

	#isoseq3 cluster ${file} ${PROCESSEDDIR}/Cluster/clustered_${sampleName}.bam --verbose --use-qvs --log-file ${PROCESSEDDIR}/Cluster/merge_${sampleName}.log

	#isoseq3 polish ${PROCESSEDDIR}/Cluster/clustered_${sampleName}.bam ${DATADIR}/${sampleName}*.subreadset.xml ${PROCESSEDDIR}/Polish/polished_${sampleName}.bam

done

## next cluster all smrt cells together
## create file with list of all filenames
ls ${PROCESSEDDIR}/Refine/*.flnc.bam > ${PROCESSEDDIR}/Refine/flnc.fofn

echo "Clustering all together"

#isoseq3 cluster ${PROCESSEDDIR}/Refine/flnc.fofn ${PROCESSEDDIR}/Cluster/clustered.bam --verbose --use-qvs --log-file ${PROCESSEDDIR}/Cluster/merge.log


#XMLFILES=($(ls ${DATADIR}/*subreadset.xml))
## to create merged xml need all the orginal files but these can be empty
## also needs .sts.xml files and merges the content of these for now copied same file across samples
#for file in ${flnc[@]}
#do
	#sample=$(basename ${file})
	#sampleName=${sample%.flnc.bam}

	#touch ${sampleName}.scraps.bam

#done

#dataset create --force --type SubreadSet --novalidate ${DATADIR}/merged.subreadset.xml ${XMLFILES[@]}

isoseq3 polish ${PROCESSEDDIR}/Cluster/clustered.bam merged.subreadset.xml ${PROCESSEDDIR}/Polish/polished.bam


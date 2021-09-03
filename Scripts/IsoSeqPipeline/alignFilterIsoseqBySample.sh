#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the high memory queue
#SBATCH --time=72:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/AlignFilterBySample.o
#SBATCH --error=LogFiles/AlignFilterBySample.e
#SBATCH --job-name=AlignFilterBySample

## this script processes each sample individually
## takes polished output
## follows pipeline at https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step but omits step to further filter isoforms based on number of full length reads

cd $SLURM_SUBMIT_DIR


source ./Config/config.txt


echo "Changing Folder to: "
echo $DATADIR

## align with minimap2
module load minimap2
fqfiles=($(ls ${PROCESSEDDIR}/Polish/*.hq.fastq.gz))
for file in ${fqfiles[@]}
do
	sample=$(basename ${file})
	sampleName=${sample#polished_}
	sampleName=${sampleName%.hq.fastq.gz}
	
	echo "Aligning " ${sampleName}

	minimap2 -ax splice -t 30 -uf --secondary=no -C5 ${REFGENOME} ${file} > ${ALIGNEDDIR}/${sampleName}.hq_isoforms.sam

	sort -k 3,3 -k 4,4n ${ALIGNEDDIR}/${sampleName}.hq_isoforms.sam > ${ALIGNEDDIR}/${sampleName}.hq_isoforms.sorted.sam
	
done


module purge
module load Miniconda2
source activate anaCogent

## collapse redundant isoforms
## to facilitate chaining output into separate folders
## need to unzip fq.gz files
if ls ${PROCESSEDDIR}/Polish/*.hq.f*q.gz 1> /dev/null 2>&1; 
then
	gunzip ${PROCESSEDDIR}/Polish/*.hq.f*q.gz
fi

samfiles=($(ls ${ALIGNEDDIR}/*.hq_isoforms.sorted.sam))
mkdir -p ${ALIGNEDDIR}/Collapsed/
for file in ${samfiles[@]}
do
	sample=$(basename ${file})
	sampleName=${sample%.hq_isoforms.sorted.sam}
	
	mkdir -p ${ALIGNEDDIR}/Collapsed/${sampleName}/

	collapse_isoforms_by_sam.py --input ${PROCESSEDDIR}/Polish/*${sampleName}*.hq*q --fq \
	   -s ${file} --dun-merge-5-shorter -o ${ALIGNEDDIR}/Collapsed/${sampleName}/out
	   
	## generate counts
	get_abundance_post_collapse.py ${ALIGNEDDIR}/Collapsed/${sampleName}/out.collapsed ${PROCESSEDDIR}/Cluster/clustered_${sampleName}.cluster_report.csv

	## Filter away 5' degraded isoforms
	filter_away_subset.py ${ALIGNEDDIR}/Collapsed/${sampleName}/out.collapsed
done


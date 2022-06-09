## bash script to automate downstream filtering of aligned isoseq transcripts
## assumes cupcake has been installed/loaded
## these steps are designed to be run on each SMRT cell separately
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires .hq.fastq.gz file located in the ${PROCESSEDDIR}/Polish folder and .hq_isoforms.sorted.sam in the ALIGNEDDIR
## follows pipeline at https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step but omits step to further filter isoforms based on number of full length reads

p=$1
basename=${p%.subreads.bam}
echo "Post alignment filtering " ${basename}

## create output folder for this sample
mkdir -p ${ALIGNEDDIR}/Collapsed/${basename}/

## may need to unzip fq.gz files
if [ -f ${PROCESSEDDIR}/Polish/polished_${basename}.hq.f*q.gz ] ## if final output file doesn't exist, run it through this loop
  then
  gunzip ${PROCESSEDDIR}/Polish/polished_${basename}.hq.f*q.gz
fi


## collapse redundant isoforms
collapse_isoforms_by_sam.py --input ${PROCESSEDDIR}/Polish/polished_${basename}.hq*q --fq -c 0.85 \
    -s ${ALIGNEDDIR}/${basename}.hq_isoforms.sorted.sam --dun-merge-5-shorter -o ${ALIGNEDDIR}/Collapsed/${basename}/out
	   
## generate counts
get_abundance_post_collapse.py ${ALIGNEDDIR}/Collapsed/${basename}/out.collapsed ${PROCESSEDDIR}/Cluster/clustered_${basename}.cluster_report.csv

## Filter away 5' degraded isoforms
filter_away_subset.py ${ALIGNEDDIR}/Collapsed/${basename}/out.collapsed


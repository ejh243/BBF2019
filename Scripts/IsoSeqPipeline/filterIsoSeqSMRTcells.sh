## bash script to automate downstream filtering of aligned isoseq transcripts
## assumes cupcake has been installed/loaded
## these steps are designed to be run on each SMRT cell separately
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires .hq.fastq.gz file located in the ${PROCESSEDDIR}/Polish folder and .hq_isoforms.sorted.sam in the ALIGNEDDIR
## follows pipeline at https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step but omits step to further filter isoforms based on number of full length reads

p=$1
basename=${p%.subreads.bam}
echo "Aligning " ${basename}

## create output folder for this sample
mkdir -p ${ALIGNEDDIR}/Collapsed/${basename}/

## need to unzip fq.gz files
gunzip ${PROCESSEDDIR}/Polish/${basename}.hq.f*q.gz

## collapse redundant isoforms
collapse_isoforms_by_sam.py --input ${PROCESSEDDIR}/Polish/*${basename}*.hq*q --fq \
    -s ${ALIGNEDDIR}/${basename}.hq_isoforms.sorted.sam --dun-merge-5-shorter -o ${ALIGNEDDIR}/Collapsed/${basename}/out
	   
## generate counts
get_abundance_post_collapse.py ${ALIGNEDDIR}/Collapsed/${basename}/out.collapsed ${PROCESSEDDIR}/Cluster/clustered_${basename}.cluster_report.csv

## Filter away 5' degraded isoforms
filter_away_subset.py ${ALIGNEDDIR}/Collapsed/${basename}/out.collapsed


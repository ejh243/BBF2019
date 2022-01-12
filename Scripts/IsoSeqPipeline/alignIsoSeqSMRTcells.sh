## bash script to automate alignment of isoseq data
## assumes minimap2 has been installed/loaded
## these steps are designed to be run on each SMRT cell separately
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires .hq.fastq.gz file are located in the ${PROCESSEDDIR}/Polish folder


p=$1
basename=${p%.subreads.bam}
echo "Aligning " ${basename}

minimap2 -ax splice -t 30 -uf --secondary=no -C5 ${REFGENOME} ${PROCESSEDDIR}/Polish/${basename}.hq.fastq.gz > ${ALIGNEDDIR}/${basename}.hq_isoforms.sam

sort -k 3,3 -k 4,4n ${ALIGNEDDIR}/${basename}.hq_isoforms.sam > ${ALIGNEDDIR}/${basename}.hq_isoforms.sorted.sam
	

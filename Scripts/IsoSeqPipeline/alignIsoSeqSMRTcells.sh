## bash script to automate processing of isoseq3 pipeline
## assumes isoseq3 has been installed/loaded
## these steps need to be run on each file separately
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires ..hq.fastq.gz files are located in the ${PROCESSEDDIR}/Polish folder
## follows pipeline at https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake:-supporting-scripts-for-Iso-Seq-after-clustering-step but omits step to further filter isoforms based on number of full length reads

p=$1
basename=${p%.subreads.bam}
echo "Aligning " ${basename}

minimap2 -ax splice -t 30 -uf --secondary=no -C5 ${REFGENOME} ${PROCESSEDDIR}/Polish/${basename}.hq.fastq.gz > ${ALIGNEDDIR}/${basename}.hq_isoforms.sam

sort -k 3,3 -k 4,4n ${ALIGNEDDIR}/${basename}.hq_isoforms.sam > ${ALIGNEDDIR}/${basename}.hq_isoforms.sorted.sam
	

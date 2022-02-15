## bash script to automate processing of isoseq3 pipeline
## assumes isoseq3 has been installed/loaded
## these steps need to be run on each file separately
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires .subread.bam, .subreads.bam.pbi, and .subreadset.xml files are located in the DATADIR

cd $DATADIR/

p=$1

echo "Processing " ${p}

basename=${p%.subreads.bam}
if [ ! -f ${PROCESSEDDIR}/CCS/${basename}.ccs.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - Circular Consensus Sequence calling"
  ## Circular Consensus Sequence calling
  ccs ${p} ${PROCESSEDDIR}/CCS/${basename}.ccs.bam --min-rq 0.9 --minPasses 1 --reportFile ${PROCESSEDDIR}/CCS/${basename}.ccs_report.txt
  
else
		echo "File Found - skipping Circular Consensus Sequence calling"
fi

if [ ! -f ${PROCESSEDDIR}/Lima/${basename}.fl.*_5p--NEB_Clontech_3p.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - Primer removal and demultiplexing"
 	
  ## Primer removal and demultiplexing
  lima --isoseq --peek-guess --dump-clips -j 24 ${PROCESSEDDIR}/CCS/${basename}.ccs.bam ${RESOURCESDIR}/primer.fasta ${PROCESSEDDIR}/Lima/${basename}.fl.bam 
  
else
		echo "File Found - skipping primer removal and demultiplexing"
fi


if [ ! -f ${PROCESSEDDIR}/Refine/${basename}.flnc.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - refine"
 
  ## refine
  isoseq3 refine --require-polya ${PROCESSEDDIR}/Lima/${basename}.fl.*_5p--NEB_Clontech_3p.bam ${RESOURCESDIR}/primer.fasta ${PROCESSEDDIR}/Refine/${basename}.flnc.bam
  
else
		echo "File Found - skipping refine step"
fi

if [ ! -f ${PROCESSEDDIR}/Cluster/clustered_${basename}.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - Clustering "

	isoseq3 cluster ${PROCESSEDDIR}/Refine/${basename}.flnc.bam ${PROCESSEDDIR}/Cluster/clustered_${basename}.bam --verbose --use-qvs --singletons --log-file ${PROCESSEDDIR}/Cluster/clustered_${basename}.log
  
else
		echo "File Found - skipping clustering"
fi

if [ ! -f ${PROCESSEDDIR}/Polish/polished_${basename}.singletons.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - Polishing " 

	isoseq3 polish ${PROCESSEDDIR}/Cluster/clustered_${basename}.bam ${DATADIR}/${basename}*.subreadset.xml ${PROCESSEDDIR}/Polish/polished_${basename}.bam
	
	isoseq3 polish ${PROCESSEDDIR}/Cluster/clustered_${basename}.singletons.bam ${DATADIR}/${basename}*.subreadset.xml ${PROCESSEDDIR}/Polish/polished_${basename}.singletons.bam

else
	echo "File Found - skipping polishing"
fi



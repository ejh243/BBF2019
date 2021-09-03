## bash script to automate processing of isoseq3 pipeline
## assumes isoseq3 has been installed/loaded
## these steps need to be run on each file separately
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires .subread.bam, .subreads.bam.pbi, and .subreadset.xml files are located in the DATADIR

cd $DATADIR/

p=$1

echo "Processing " ${p}

basename=${p%.subreads.bam}
if [ ! -f ${PROCESSEDDIR}/Refine/${basename}.flnc.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - processing"
  ## Circular Consensus Sequence calling
  ccs ${p} ${PROCESSEDDIR}/CCS/${basename}.ccs.bam --min-rq 0.9 --minPasses 1 --reportFile ${PROCESSEDDIR}/CCS/${basename}.ccs_report.txt
	
  ## Primer removal and demultiplexing
  lima --isoseq --peek-guess --dump-clips -j 24 ${PROCESSEDDIR}/CCS/${basename}.ccs.bam ${RESOURCESDIR}/primer.fasta ${PROCESSEDDIR}/Lima/${basename}.fl.bam 
		  
  ## refine
  isoseq3 refine --require-polya ${PROCESSEDDIR}/Lima/${basename}.fl.Clontech_5p--NEB_Clontech_3p.bam ${RESOURCESDIR}/primer.fasta ${PROCESSEDDIR}/Refine/${basename}.flnc.bam
  
else
		echo "File Found - skipping"
fi






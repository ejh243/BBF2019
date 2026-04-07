## bash script to automate processing of isoseq3 pipeline 
## Updated to use the latest Isoseq3 tools - version 4.3.0 (https://isoseq.how/getting-started.html)
## assumes isoseq3 has been installed/loaded 
## these steps need to be run on each file separately
## do not store any sensitive data use config file to specify filepaths etc. 
## this script requires .subread.bam, .subreads.bam.pbi, and .subreadset.xml files are located in the DATADIR

p=$1

echo "Processing " ${p}
basename=$(basename "$p" .subreads.bam)
echo "Basename: $basename"

echo "Processing ${p}"
echo "PWD: $(pwd)"
echo "DATADIR: ${DATADIR}"
echo "Input file: ${p}"
echo "Basename: ${basename}"
ls -lh "${p}"
echo "Batch script DATADIR: ${DATADIR}"

# Step 1: CCS - Generate circular consensus sequences (ccs) from subreads
if [ ! -f ${PROCESSEDDIR}/CCS/${basename}.ccs.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - Circular Consensus Sequence calling"
  ## Circular Consensus Sequence calling
  ccs ${p} ${PROCESSEDDIR}/CCS/${basename}.ccs.bam --min-rq 0.9 --min-passes 1 --report-file ${PROCESSEDDIR}/CCS/${basename}.ccs_report.txt
  
else
		echo "File Found - skipping Circular Consensus Sequence calling"
fi


# Step 2: Lima - Remove cDNA primers and demultiplexing barcoded data 
if [ ! -f ${PROCESSEDDIR}/Lima/${basename}.fl.*_5p--NEB_Clontech_3p.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - Primer removal and demultiplexing"
 	
  ## Primer removal and demultiplexing
  lima --isoseq --peek-guess --dump-clips --num-threads 24 ${PROCESSEDDIR}/CCS/${basename}.ccs.bam ${RESOURCESDIR}/primer.fasta ${PROCESSEDDIR}/Lima/${basename}.fl.bam 
  
else
		echo "File Found - skipping primer removal and demultiplexing"
fi


# Step 3: Refine - Remove polyA and concatemers from FL reads and generate FLNC transcripts
if [ ! -f ${PROCESSEDDIR}/Refine/${basename}.flnc.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - refine"
 
  ## refine
  isoseq refine --require-polya ${PROCESSEDDIR}/Lima/${basename}.fl.*_5p--NEB_Clontech_3p.bam ${RESOURCESDIR}/primer.fasta ${PROCESSEDDIR}/Refine/${basename}.flnc.bam
  
else
		echo "File Found - skipping refine step"
fi


# Step 4: Cluster2 - Cluster FLNC reads and generate transcripts
# Note - polish step not required in newer pipeline / cluster2 tool 
if [ ! -f ${PROCESSEDDIR}/Cluster2/clustered_${basename}.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "File not found - Clustering "

	isoseq cluster2 ${PROCESSEDDIR}/Refine/${basename}.flnc.bam ${PROCESSEDDIR}/Cluster2/clustered_${basename}.bam --singletons --log-file ${PROCESSEDDIR}/Cluster/clustered_${basename}.log
  
else
		echo "File Found - skipping clustering"
fi

# End of bulk Iso-Seq workflow. Next to continue to pigeon workflow. 


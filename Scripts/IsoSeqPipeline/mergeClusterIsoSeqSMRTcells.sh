## bash script to automate merging and clustering of processed IsoSeq data
## Updated to use the latest Isoseq3 tools - version 4.3.0 (https://isoseq.how/getting-started.html)
## assumes isoseq3 has been installed/loaded 
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires .flnc.bam files are located in the ${PROCESSEDDIR}/Refine folder

# If multiple SMRTcells found in PROCESSEDDIR, then create a list of all <movie>.flnc.bam using
# This list can be used as input for cluster2 step


p=$1

basename=$(basename "$p" .subreads.bam)

echo "Looking for multiple SMRTcells in: ${PROCESSEDDIR}"

if [ ! -f ${PROCESSEDDIR}/Refine/flnc.fofn ] ## if final output file doesn't exist, run it through this loop
  then
  echo "Ouput file not found - Running Merging step"

	ls ${PROCESSEDDIR}/Refine*.flnc.bam > flnc.fofn
  
else
		echo "Ouput file Found - skipping merging"
fi
 

# Cluster2 - Cluster FLNC reads and generate transcripts
# Note - polish step not required in newer pipeline / cluster2 tool 

if [ ! -f ${PROCESSEDDIR}/Cluster2/clustered_${basename}.bam  ## if final output file doesn't exist, run it through this loop
  then
  echo "Ouput file not found - Running Clustering step"

	isoseq cluster2 ${PROCESSEDDIR}/Refine/${basename}.flnc.bam ${PROCESSEDDIR}/Cluster2/clustered_${basename}.bam --singletons --log-file ${PROCESSEDDIR}/Cluster/clustered_${basename}.log
  
else
		echo "Ouput file Found - skipping clustering"
fi
 

# End of script 
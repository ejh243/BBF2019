## bash script to automate processing of isoseq3 pipeline 
## Updated to use the latest Isoseq3 tools - version 4.3.0 (https://isoseq.how/getting-started.html)
## assumes isoseq3 has been installed/loaded 
## these steps need to be run on each file separately
## do not store any sensitive data use config file to specify filepaths etc. 
## this script requires .subread.bam, .subreads.bam.pbi, and .subreadset.xml files are located in the DATADIR

p=$1

echo "Processing " ${p}
basename=$(basename "$p" .subreads.bam)

echo "DATADIR: ${DATADIR}"
echo "Input file: ${p}"

# Step 1: CCS - Generate circular consensus sequences (ccs) from subreads
if [ ! -f ${PROCESSEDDIR}/CCS/${basename}.ccs.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "Ouput file not found - Running Circular Consensus Sequence calling"
  # Run Circular Consensus Sequence calling
  ccs ${p} ${PROCESSEDDIR}/CCS/${basename}.ccs.bam --min-rq 0.9 --min-passes 1 --report-file ${PROCESSEDDIR}/CCS/${basename}.ccs_report.txt
  
else
		echo "Ouput file Found - skipping Circular Consensus Sequence calling"
fi


# Step 2: Lima - Remove cDNA primers and demultiplexing barcoded data 
if [ ! -f ${PROCESSEDDIR}/Lima/${basename}.fl.*_5p--NEB_Clontech_3p.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "Ouput file not found - Running Primer removal and demultiplexing"
  # Run Primer removal and demultiplexing
  lima --isoseq --peek-guess --dump-clips --num-threads 24 ${PROCESSEDDIR}/CCS/${basename}.ccs.bam ${RESOURCESDIR}/primer.fasta ${PROCESSEDDIR}/Lima/${basename}.fl.bam 
  
else
		echo "Ouput file Found - skipping primer removal and demultiplexing"
fi


# Step 3: Refine - Remove polyA and concatemers from FL reads and generate FLNC transcripts
if [ ! -f ${PROCESSEDDIR}/Refine/${basename}.flnc.bam ] ## if final output file doesn't exist, run it through this loop
  then
  echo "Ouput file not found - Running Isoseq refine"
  # Run Refine
  isoseq refine --require-polya ${PROCESSEDDIR}/Lima/${basename}.fl.*_5p--NEB_Clontech_3p.bam ${RESOURCESDIR}/primer.fasta ${PROCESSEDDIR}/Refine/${basename}.flnc.bam
  
else
		echo "Ouput file Found - skipping refine step"
fi

echo "End of Isoseq per SMRT cell bulk processing script."
# End of script 
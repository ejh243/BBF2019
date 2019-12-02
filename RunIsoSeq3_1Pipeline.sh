## bash script to automate processing of isoseq3.1 pipeline
## assumes isoseq3.1 has been installed
## do not store any sensitive data use config file to specify filepaths etc.
## this script requires .subread.bam, .subreads.bam.pbi, and .subreadset.xml files are located in the DATADIR

cd $DATADIR/


## run first steps on each sample individually
mkdir -p CCS
mkdir -p Lima
mkdir -p Refine

## this command can be used to process all relevant files in DATADIR 
## if only a subset are created provide list in FilesToProcess.txt and hash out line below.
cd RawData/
ls *.subreads.bam > FilesToProcess.txt
cd ..

## output version of ccs
ccs --version
## output version of lima
lima --version

while read p; do
  basename=${p%.subreads.bam}
	if [ -f Refine/${basename}.flnc.bam ] ## if final output file doesn't exist, run it through this loop
	  then
		else    
	  ## Circular Consensus Sequence calling
	  ccs RawData/${p} CCS/${basename}.ccs.bam --min-rq 0.9 --minPasses 1 --reportFile CCS/${basename}.ccs_report.txt
		
	  ## Primer removal and demultiplexing
	  lima CCS/${basename}.ccs.bam Resources/primer.fasta Lima/${basename}.fl.bam --isoseq --no-pbi --peek-guess --dump-clips -j 24
		
	  ## refine
	  isoseq3 refine --require-polya Lima/${basename}.fl.Clontech_5p--NEB_Clontech_3p.bam Resources/primer.fasta Refine/${basename}.flnc.bam
	fi
done < RawData/FilesToProcess.txt



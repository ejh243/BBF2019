## merge multiple data across samples and continue processing
cd $DATADIR/

allBam=$(ls Refine/*.flnc.bam)

mkdir -p Merged
dataset create --type TranscriptSet Merged/merged.flnc.xml ${allBam}


## cluster 
isoseq3 cluster Merged/merged.flnc.xml Cluster/polished.bam --verbose --use-qvs
isoseq3 summarize Cluster/polished.bam summary.csv



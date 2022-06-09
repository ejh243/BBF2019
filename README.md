# Repository of stripts for processing PacBio IsoSeq data


## Data structure

Datasets are organised by processing stage and then into folders by project name. 

```
project
│   README.md
│   file001.txt    
│
└───RawData
   │
   └───PacBio

   └───RNASeq
       └───Project1
           │   sample1_R1.fastq.gz
           │   sample1_R2.fastq.gz
           │   ...       
           └───fastqc
               | sample1_fastqc.zip
               | sample1_fastqc.html
               | ...
               | multiqc_report.html
               └───multiqc_data 
           └───trimmed
               | sample1_R1_val_1.fq.gz
               | sample1_R2_val_2.fq.gz
               | sample1_R1_val_1_fastqc.zip
               | sample1_R2_val_2_fastqc.zip 
               | sample1_R1_val_1_fastqc.html
               | sample1_R2_val_2_fastqc.html            
               | sample1_R1_val_1.fq.gz_trimming_report.txt
               | sample1_R1_val_1.fq.gz_trimming_report.txt               
               | ...          
   └───Aligned
    │
    └───RNASeq
       └───Project1
           | sample1_Aligned.sortedByCoord.out.bam
           | sample1_Log.final.out
           | sample1_Log.progress.out
           | sample1_Log.out
           | sample1_SJ.out.tab           
           | ...       
           | multiqc_report.html
           └───rnaseqc             
           └───multiqc_data           
                     
└───GeneCounts
    │   file021.txt
    │   file022.txt 
└───Scripts
    │   file021.txt
    │   file022.txt    
```


Details for running the isoseq pipeline



Run the following commands from the github repository main folder

## SMRT cell processing into master transcriptome

1. sbatch Scripts/JobSubmission/batchProcessIsoSeqData.sh
	searchs for all files with sufix .subreads.bam then executes the following scripts for each sample:
	
	* processIsoSeqSMRTcells.sh runs the isoseq3 pipeline including ccs, lima, refine, cluster and polish steps.
	* alignIsoSeqSMRTcells.sh aligns the processed data using minimap2
	* filterIsoSeqSMRTcells.sh performs filtering on the on the aligned transcripts using cupcake scripts
	* alignShortReadSMRTcells.sh aligns short read data to gff for use in sqanti qc
	* sqanti3SMRTcells.sh QCs SMRT-cell level transcriptomes prior to merging
	
2. smrt cell level qc summary

3. sbatch Scripts/IsoSeqPipeline/createMasterTranscriptomeFromSMRTcells.sh
	use talon to mark reads for internal priming and merge smrt cell transcriptomes into single master transcriptomes. Compare to GENCODE for annotation of known isoforms

4. sbatch Scripts/IsoSeqPipeline/.sh
	use talon to annotate master transcriptome with 

4. characterise transcriptome

## RNA-Seq data pipeline

1. sbatch --array=<number of jobs> Scripts/JobSubmission/runRNASeqAlignment.sh <folder with fastq files> <Project Name>
	* Scripts/RNASeq/alignShortReadData.sh: fastqc, trimGalore, align with STAR to GENCODE
	* Scripts/RNASeq/rsemGENCODE.sh: quantify GENCODE transcripts with RSEM
	* Scripts/RNASeq/rnaseqQC.sh: RNASEQ-QC
	
2. sbatch --array=0-49%10 Scripts/JobSubmission/runAlignmentBrainTranscriptome.sh <folder with fastq files> <Project Name>
	* Scripts/RNASeq/rsemBrainTranscriptome.sh: quantify Brain Transcriptome with RSEM
	
3. sbatch JobSubmission/featureCountsGENCODE.sh <Project Name>
	* uses featureCounts to calculate exon and gene level counts for GENCODE transcriptome
	
	
## Differential Expression analyses

1. Rscript Scripts/DiffSplicing/edgeRExpression.r <directory with RSEM counts> <sample sheet> <output directory>


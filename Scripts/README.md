# Iso-Seq updated pipeline 
Runs the latest Isoseq3 tools - version 4.3.0 (https://isoseq.how/getting-started.html)

The required software tools are loaded on a conda environment 
created using Scripts/Utilities/createCondaEnv_updated.sh
load using `source activate isoseq_tools` 


## Details for running the isoseq pipeline 
1. Per-SMRT cell preprocessing (ccs, lima and refine) <br>
    This creates HiFi reads, removes primers and refine into full length non-concatamer reads <br>
    Designed to run on each SMRT cell individually (parallelisation/batch processing) <br>
    `sbatch Scripts/IsoSeqPipeline/01_preprocess_IsoSeqSMRTcells.sh` <br>

2. Compilation of SMRT cell processing QC <br>
    This parses relevant QC metrics from CCS and Refine reports and compiles into a TSV
    `sbatch Scripts/IsoSeqPipeline/02_compile_SMRTcell_QCreports.sh` <br>

3. Preparation of input files before merging data across SMRT cells 
    Includes editing BAM headers (SM tag) to backtrace each SMRT run after merging <br>
    `sbatch Scripts/IsoSeqPipeline/03_prepare_merged_SMRTcell_input.sh`

4.  Build master transcriptome using data across all SMRT cells (cluster2, pbmm2, collapse) <br>
    Clusters similar transcripts, aligns clustered FLNC reads to reference genome and <br>
    collapses FLNC transcripts into unique isoforms based on genomic locus and exonic structure similarity <br>
    `sbatch Scripts/IsoSeqPipeline/04_build_master_transcriptome.sh`

5.  Classify and filter FLNC transcript isoforms into categories against a reference annotation <br>
    Uses Pigeon - PacBio's new Transcript Toolkit (https://isoseq.how/classification/pigeon.html) <br>
    Isoform categories (Based off SQANTI3 - FSM, ISM, NIC/NNC, genetic) <br>
    `sbatch Scripts/IsoSeqPipeline/05_annotate_master_transcriptome.sh`



# RNA-Seq pipeline 
Runs standard RNA-Seq pipeline for alignment and quantification of gene expression using short reads <br>
load using `source activate rnaseq_tools` <br>

## Scripts used so far:
1. Preprocessing of short reads (fastqc, trim_galore) <br>
    `sbatch Scripts/RNASeq/01_preprocess_ShortReads.sh`<br>

2. Alignment against reference with GENCODE (STAR) <br>
    `sbatch Scripts/RNASeq/02_align_ShortReads.sh` <br>

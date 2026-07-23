# Iso-Seq and RNA-seq Analysis Pipelines
This repository contains analysis pipelines for processing PacBio Iso-Seq and RNA-seq data. <br> 

The Iso-Seq pipeline implements the current PacBio IsoSeq3 workflow (v4.3.0) to generate a high-confidence, non-redundant transcriptome from one or more SMRT cells. The RNA-seq pipeline is currently under development and will provide complementary gene and transcript expression analyses<br> <br>

### Repository structure
```
.
├── Config/
│   └── config.txt
├── Scripts/
│   ├── IsoSeqPipeline/
│   ├── RNASeq/
│   ├── Utilities/
│   ├── JobSubmission/
│   └── tamaPipeline/
└── README.md
```

**Config:** The `Config/` directory contains the project configuration file, which can be used to specify data locations as well as input/output directories. This allows the pipeline scripts to remain unchanged between projects. 

**Software Installation:** 
Required software is installed using the scripts in `Scripts/Utilities/` 
Running `Scripts/Utilities/createCondaEnv_updated.sh` creates two Conda environments: isoseq_tools and rnaseq_tools. The required environment is activated automatically within the pipeline scripts, so manual activation is generally not required. 
<br> The Conda environment (.yml) files are also located here <br> <br>

## Iso-Seq pipeline 
An automated PacBio Iso-Seq analysis pipeline for constructing and annotating a master transcriptome from multiple SMRT cells <br>
### Details for running the isoseq pipeline 
**1. Per-SMRT cell preprocessing** <br>
    Runs CCS, Lima and Refine on each SMRT cell independently.
    This step generates HiFi (CCS) reads, removes primer sequences and produces full-length non-concatemer (FLNC) reads <br> <br>
    Designed to be run separately for each SMRT cell, allowing parallel or batch processing <br>
    `sbatch Scripts/IsoSeqPipeline/01_preprocess_IsoSeqSMRTcells.sh` <br>

**2. SMRT cell QC report compilation** <br>
    Parses QC metrics from the CCS and Refine reports and compiles them into a single tab-separated summary. <br> Runs `parse_SMRTcell_QC_metrics.py` through below bash wrapper script. <br> 
    `sbatch Scripts/IsoSeqPipeline/02_compile_SMRTcell_QCreports.sh` <br>

**3. Data preparation prior to merging across SMRT cells** <br>
    This step edits BAM headers by adding or updating the `SM` tag, enabling reads to be traced back to their original SMRT cell after merging. <br>
    `sbatch Scripts/IsoSeqPipeline/03_prepare_merged_SMRTcell_input.sh`

**4.  Build the master transcriptome** <br>
    Processes data across all SMRT cells to construct a unified transcriptome. <br> This step: clusters similar FLNC reads (cluster2)
    aligns clustered reads to the reference genome (pbmm2)
    collapses redundant transcripts into unique isoforms (collapse) <br> 
    Isoforms are merged based on genomic location and exon structure. <br>
    `sbatch Scripts/IsoSeqPipeline/04_build_master_transcriptome.sh`

**5.  Annotate the master transcriptome** <br> 
    Classifies transcript isoforms relative to a reference annotation using PacBio's Pigeon Transcript toolkit  <br>
    Assigns transcript categories based on the SQANTI3 classification framework, (FSM, ISM, NIC/NNC, genic and other). More information available at: https://isoseq.how/classification/pigeon.html <br>
    `sbatch Scripts/IsoSeqPipeline/05_annotate_master_transcriptome.sh`

### Example project layout
A typical project will generate an output directory structure similar to: 
```
Project/ 
├── RawData/ 
│   └── PacBio/ 
├── 01_IsoSeqPreprocessed/
│   ├── CCS/
│   ├── Lima/
│   ├── Refine/
│   └── QC
├── 02_IsoSeqMerged/
│   ├── PreparedFLNC/
│   ├── Clustered/
│   └── Aligned/
├── 03_MasterTranscriptome/
│   ├── Isoforms/
│   └── Annotation/
```

### Quick start
1. Clone this repository.
2. Create the Conda environments.
3. Edit `Config/config.txt` for your project.
4. Submit the pipeline scripts in numerical order.


## RNA-Seq pipeline 
Runs standard RNA-Seq pipeline for alignment and quantification of gene expression using short reads <br>
load using `source activate rnaseq_tools` <br>

### Scripts used so far:
1. Preprocessing of short reads (fastqc, trim_galore) <br>
    `sbatch Scripts/RNASeq/01_preprocess_ShortReads.sh`<br>

2. Alignment against reference with GENCODE (STAR) <br>
    `sbatch Scripts/RNASeq/02_align_ShortReads.sh` <br>

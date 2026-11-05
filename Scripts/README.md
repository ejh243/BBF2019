# Iso-Seq and RNA-seq Analysis Pipelines
This repository contains analysis pipelines for processing PacBio Iso-Seq and RNA-seq data. <br> 

The **Iso-Seq pipeline** implements the current PacBio IsoSeq3 workflow (v4.3.0) to generate a high-confidence, non-redundant transcriptome from one or more SMRT cells. The **RNA-seq pipeline** is currently under development and will provide complementary gene and transcript expression analyses<br> <br>

### Repository structure
```
.
├── Config/
│   └── config.txt
├── Scripts/
│   ├── IsoSeqPipeline/
│   ├── IsoSeqAnalysis/
│   ├── RNASeqPipeline/
│   ├── RNASeqAnalysis/
│   └── Utilities/
└── README.md
```

**Config:** The `Config/` directory contains the project configuration file, which can be used to specify data locations as well as input/output directories. This allows the pipeline scripts to remain unchanged between projects. 

**Software Installation:** 
Required software is installed using the scripts in `Scripts/Utilities/` 
Running the `createCondaEnv.sh` script creates two custom Conda environments: isoseq_tools and rnaseq_tools. The script also contains instructions to install Sqanti on a separate conda environment. <br> The required environment is activated automatically within the pipeline scripts, so manual activation is generally not required. 
<br> The Conda environment (.yml) files are also located within this directory. <br> <br>

## 1. Iso-Seq pipeline 
An automated PacBio Iso-Seq analysis pipeline for constructing and annotating a master transcriptome from multiple SMRT cells <br>
### [A] Details for running the isoseq pipeline 
**1. Per-SMRT cell preprocessing** <br>
    Runs on each SMRT cell independently using SLURM array for batch processing. <br> This script: <br>
    - generates HiFi reads from raw subreads (css) <br> 
    - demultiplexes and removes primer sequences (lima) <br> 
    - produces full-length non-concatemer "FLNC" reads (isoseq refine) <br> <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/01_preprocess_IsoSeqSMRTcells.sh` <br> 

**2. SMRT cell QC report compilation** <br>
    Parses QC metrics from the CCS and Refine reports and compiles them into a single tab-separated summary. Acts as a wrapper script to run custom script `parse_SMRTcell_QC_metrics.py` <br> <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/02_compile_SMRTcell_QCreports.sh` <br> 

**3. Data preparation prior to merging across SMRT cells** <br>
    Edits BAM headers by adding and updating the `SM` tag, enabling reads to be traced back to their original SMRT cell after merging. Acts as a preliminary step to enable merging of FLNC data across SMRT cells. <br> <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/03_prepare_merged_SMRTcell_input.sh`

**4.  Build the master transcriptome** <br>
    Processes data across all SMRT cells to construct a unified transcriptome. <br> This script: <br>
    - clusters sequence-similar FLNC transcripts (isoseq cluster2) <br>
    - maps clustered transcripts to a reference genome (pbmm2) <br>
    - collapses redundant transcripts into unique isoforms based on exonic structure and genomic coordinates (isoseq collapse) <br> <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/04_build_master_transcriptome.sh`

**5.  Annotate the master transcriptome** <br> 
    The master transcriptome is annotated using a combined SQANTI3 and PacBio Pigeon workflow. This approach uses SQANTI to incorporate short-read junction coverage information from RNA-seq data and uses Pigeon workflow for cleaner filtering and QC. <br> *This is a workaround for a current issue with Pigeon's handling of short-read junction coverage information, as confirmed by PacBio.* <br> 

**5a. SQANTI3 isoform classification:** 
    Classifies collapsed transcript isoforms relative to a reference annotation, using the SQANTI structural classification framework (FSM, ISM, NIC, NNC and other categories). More information available at: https://isoseq.how/classification/pigeon.html <br> <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/05_classify_master_transcriptome.sh` <br> 

**5b. Pigeon classification, coverage patching and filtering:** <br>
    This script: <br>
    - Runs PacBio Pigeon independently to generate its transcript classification and associated junctions files <br> 
    - Updates the resulting classification with the junction coverage metrics obtained from SQANTI3, using the custom script `patch_coverage_data_into_pigeon.py`. *This is due to Pigeon being unable to directly parse the SQANTI3 classification file* <br>
    - Filters isoforms from artifacts using Pigeon's modernised filtering functionality <br> <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/05_annotate_master_transcriptome.sh`

### [B] Example project layout
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
│   └── Mapped/
├── 03_MasterTranscriptome/
│   ├── Isoforms/
│   └── Annotation/
        ├── SQANTI_classification/
        ├── Pigeon_classification/
        └── Final_isoform_classification/
```

## 2. RNA-Seq pipeline 
Runs standard RNA-Seq pipeline for alignment and quantification of gene expression using short reads <br>
load using `source activate rnaseq_tools` <br>

### Scripts used so far:
1. Preprocessing of short reads (fastqc, trim_galore) <br>
    `sbatch Scripts/RNASeq/01_preprocess_ShortReads.sh`<br>

2. Alignment against reference with GENCODE (STAR) <br>
    `sbatch Scripts/RNASeq/02_align_ShortReads.sh` <br>


### Quick start
1. Clone this repository.
2. Create the Conda environments.
3. Edit `Config/config.txt` for to specify project directories, reference files and other required input paths.
4. Submit the pipeline scripts in numerical order.
The pipeline scripts are designed for SLURM submission and should be run from the main repository directory.

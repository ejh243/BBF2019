# Iso-Seq and RNA-seq Analysis Pipelines
This repository contains analysis pipelines for processing PacBio Iso-Seq and RNA-seq data. <br> 

The **Iso-Seq pipeline** implements the current PacBio IsoSeq3 workflow (v4.3.0) to generate a high-confidence, non-redundant transcriptome from one or more SMRT cells. The **RNA-seq pipeline** provides complementary gene and transcript quantification analyses using short read sequencing data. <br> <br>

## Repository structure
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
    - produces full-length non-concatemer "FLNC" reads (isoseq refine) <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/01_preprocess_IsoSeqSMRTcells.sh` <br> 

**2. SMRT cell QC report compilation** <br>
    Parses QC metrics from the CCS and Refine reports and compiles them into a single tab-separated summary. Acts as a wrapper script to run custom script `parse_SMRTcell_QC_metrics.py` <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/02_compile_SMRTcell_QCreports.sh` <br> 

**3. Data preparation prior to merging across SMRT cells** <br>
    Edits BAM headers by adding and updating the `SM` tag, enabling reads to be traced back to their original SMRT cell after merging. Acts as a preliminary step to enable merging of FLNC data across SMRT cells. <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/03_prepare_merged_SMRTcell_input.sh`

**4.  Build the master transcriptome** <br>
    Processes data across all SMRT cells to construct a unified transcriptome. <br> This script: <br>
    - clusters sequence-similar FLNC transcripts (isoseq cluster2) <br>
    - maps clustered transcripts to a reference genome (pbmm2) <br>
    - collapses redundant transcripts into unique isoforms based on exonic structure and genomic coordinates (isoseq collapse) <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/04_build_master_transcriptome.sh`

**5.  Annotate the master transcriptome** <br> 
    The master transcriptome is annotated using a combined SQANTI3 and PacBio Pigeon workflow. This approach uses SQANTI to incorporate short-read junction coverage information from RNA-seq data and uses Pigeon workflow for cleaner filtering and QC. <br> *This is a workaround for a current issue with Pigeon's handling of short-read junction coverage information, as confirmed by PacBio.* <br> 

**5a. SQANTI3 isoform classification:** 
    Classifies collapsed transcript isoforms relative to a reference annotation, using the SQANTI structural classification framework (FSM, ISM, NIC, NNC and other categories). More information available at: https://isoseq.how/classification/pigeon.html <br>
    Usage: `sbatch Scripts/IsoSeqPipeline/05_classify_master_transcriptome.sh` <br> 

**5b. Pigeon classification, coverage patching and filtering:** <br>
    This script: <br>
    - Runs PacBio Pigeon independently to generate its transcript classification and associated junctions files <br> 
    - Updates the resulting classification with the junction coverage metrics obtained from SQANTI3, using the custom script `patch_coverage_data_into_pigeon.py`. *This is due to Pigeon being unable to directly parse the SQANTI3 classification file* <br>
    - Filters isoforms from artifacts using Pigeon's modernised filtering functionality <br>
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
│       ├── SQANTI_classification/
│       ├── Pigeon_classification/
│       └── Final_isoform_classification/
```

## 2. RNA-Seq pipeline 
This standard RNA-Seq pipeline generates two main types of expression quantification data: <br>
**Gene-level quantification:** Raw read counts generated using STAR-aligned BAM files and featureCounts. <br>
**Transcript-level quantification:** Transcript abundance estimates generated using Salmon, including estimated counts and TPM values. <br> <br>
These outputs can be used for downstream analyses, including differential gene expression analysis, transcript-level expression analysis, and functional enrichment analysis. <br>

### [A] Details for running the isoseq pipeline 
**1. Short Read preprocessing** <br>
    This script: 
    - performs initial quality assessment of raw paired-end FASTQ files (fastqc) <br> 
    - removes adapter sequences and trims low-quality bases (trim galore) <br> 
    Usage: `sbatch Scripts/RNASeqPipeline/01_preprocess_ShortReads.sh` 

**2. Reference index generation** <br>
    Prepares the reference genome and transcriptome indices required by STAR and Salmon <br>
    Usage: `sbatch Scripts/RNASeqPipeline/02_create_alignment_index.sh` 

**3. Gene-level read alignment** <br>
    Performs splice-aware alignment of trimmed RNA-seq reads against the reference genome or transcriptome (STAR) <br> 
    Usage: `sbatch Scripts/RNASeqPipeline/03_perform_read_alignment.sh`

**4. Gene-level expression quantification** <br>
    Assigns aligned reads to annotated genomic features and generates gene-level read counts (FeatureCounts) <br>
    Usage: `sbatch Scripts/RNASeqPipeline/04_quantify_gene_features.sh` 

**5. Transcript-level expression quantification** <br>
    Quantifies transcript abundance directly from the preprocessed paired-end FASTQ files (Salmon) <br> 
    Usage: `sbatch Scripts/RNASeqPipeline/05_quantify_transcripts.sh`

### [B] Example project layout
A typical project will generate an output directory structure similar to: 
```
Project/ 
├── RawData/ 
│   └── RNASeq/ 
├── 01_processed_reads/
│   ├── fastqc_raw/
│   ├── fastqc_trimmed/
│   └── reads_trimmed/
├── 02_aligned_reads/
│   ├── alignment_against_gencode/
│   └── alignment_against_transcriptome/
├── 03_gene_quantification/
│   ├── gene_counts/
│   └── analysis/
├── 04_transcript_quantification/
│   ├── salmonResults/
│   └── analysis/
```

## Quick start
1. Clone this repository.
2. Create the Conda environments.
3. Edit `Config/config.txt` for to specify project directories, reference files and other required input paths.
4. Submit the pipeline scripts in numerical order.
The pipeline scripts are designed for SLURM submission and should be run from the main repository directory.

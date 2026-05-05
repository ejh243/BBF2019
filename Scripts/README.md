# Iso-Seq updated pipeline 
Runs the latest Isoseq3 tools - version 4.3.0 (https://isoseq.how/getting-started.html)

The required software tools are loaded on a conda environment 
created using Scripts/Utilities/createCondaEnv_updated.sh
load using `source activate isoseq_tools` 


## Details for running the isoseq pipeline 
1. Per-SMRT cell preprocessing (ccs, lima and refine)
    This creates HiFi reads, removes primers and refine into full length non-concatamer reads
    Designed to run on each SMRT cell individually (parallelisation/batch processing)
    `sbatch Scripts/IsoSeqPipeline/processIsoSeqSMRTcells_v2.sh`

*Notes: After this, need to edit each BAM file to represent a different "BioSample". Could try to write a script for this.*

2.  Merge across SMRT cells and cluster sequence similar FLNC transcripts (cluster2) 
    `sbatch Scripts/IsoSeqPipeline/mergeClusterIsoSeqSMRTcells.sh`

3.  Align clustered FL transcripts to reference genome hg38 (pbmm2) 
    `sbatch Scripts/IsoSeqPipeline/alignMergedFlncReads.sh`

4. Collapse FLNC transcripts into unique isoforms (collapse)
    Collapses based on genomic locus and exonic structure similarity
    This essentially creates a master transcriptome from all SMRT cells 
    `sbatch Scripts/IsoSeqPipeline/collapseIsoforms.sh`

*Notes: Step 1 runs jobs for each individual SMRTcell whereas the scripts 2-4 could potentially be run from a single workflow / wrapper script?* 


## Pigeon pipeline 
This is a new PacBio Transcript Toolkit used to classify and filter FLNC transcript isoforms into categories against a reference annotation. (Based off SQANTI3, https://isoseq.how/classification/pigeon.html)

1. Prepare (sort and index) reference files and input transcript files for the next step.
    Can provide a list of genome annotation GTF file, the reference file and the transcript GFF file output from isoseq collapse. 
    `sbatch Scripts/PigeonPipeline/prepareInputFiles.sh`

2. Classify isoforms based on different categories (FSM, ISM, NIC/NNC, genetic)
    `sbatch Scripts/PigeonPipeline/classifyIsoforms.sh`

*Notes: The above two scripts could be put into a single script. With file checks on the "prepare" step before proceeding.* 




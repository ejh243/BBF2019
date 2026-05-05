# Iso-Seq updated pipeline 
Runs the latest Isoseq3 tools - version 4.3.0 (https://isoseq.how/getting-started.html)

The required software tools are loaded on a conda environment 
created using Scripts/Utilities/createCondaEnv_updated.sh
load using `source activate isoseq_tools` 


## Details for running the isoseq pipeline 
1. Per-SMRT cell preprocessing (ccs, lima and refine) <br>
    This creates HiFi reads, removes primers and refine into full length non-concatamer reads <br>
    Designed to run on each SMRT cell individually (parallelisation/batch processing) <br>
    `sbatch Scripts/IsoSeqPipeline/processIsoSeqSMRTcells_v2.sh`

2. Script to edit BAM headers (SM tag) to backtrace each SMRT run after merging <br>
    `sbatch Scripts/IsoSeqPipeline/changeBAMHeaders.sh`

3.  Merge across SMRT cells and cluster sequence similar FLNC transcripts (cluster2) <br>
    `sbatch Scripts/IsoSeqPipeline/mergeClusterIsoSeqSMRTcells.sh`

4.  Align clustered FL transcripts to reference genome hg38 (pbmm2) <br>
    `sbatch Scripts/IsoSeqPipeline/alignMergedFlncReads.sh`

5. Collapse FLNC transcripts into unique isoforms (collapse) <br>
    Collapses based on genomic locus and exonic structure similarity <br>
    This essentially creates a master transcriptome from all SMRT cells <br>
    `sbatch Scripts/IsoSeqPipeline/collapseIsoforms.sh`


## Pigeon pipeline 
This is a new PacBio Transcript Toolkit used to classify and filter FLNC transcript isoforms into categories against a reference annotation. (Based off SQANTI3, https://isoseq.how/classification/pigeon.html)

1. Prepare (sort and index) reference files and input transcript files for the next step.
    Can provide a list of genome annotation GTF file, the reference file and the transcript GFF file output from isoseq collapse. <br>
    `sbatch Scripts/PigeonPipeline/prepareInputFiles.sh`

2. Classify isoforms based on different categories (FSM, ISM, NIC/NNC, genetic) <br>
    `sbatch Scripts/PigeonPipeline/classifyIsoforms.sh`


*Notes: Could merge some of the scripts or create wrapper scripts. Potentially combine steps 3-5 of isoseq pipeline and leave script 2 as optional? Could also combine the two pigeon scripts into a single script.*





#!/usr/bin/env bash
set -euo pipefail

## set up conda environment for Modern IsoSeq and RNA-Seq pipelines 

## IsoSeq environment 
#   Based on IsoSeq CLI workflow (version 4.3.0)
#   No Cupcake scripts, uses Python 3 only
#   Also installing Python data science software for downstream QC 
conda create -n isoseq_tools -y \
    -c conda-forge \
    -c bioconda \
    --strict-channel-priority \
    python=3.10 \
    isoseq \
    lima \
    pbmm2 \
    pbpigeon \
    pbbam \
    pbccs \
    samtools \
    bamtools \
    ipykernel \
    jupyter \
    pandas \
    numpy \
    matplotlib \
    seaborn 


## RNASeq environment 
#   Installing gene-level and transcript-level quantification tools 
#   And some R packages for downstream analysis 
conda create -n rnaseq_tools -y \
    -c conda-forge \
    -c bioconda \
    --strict-channel-priority \
    python=3.10 \
    fastqc \
    star \
    trim-galore \
    salmon \
    kallisto\
    samtools \
    r-essentials \
    bioconductor-deseq2 \
    subread 


# To create environment.yml (with software versions)
# conda env export --no-builds -n isoseq_tools > isoseq_tools_env_full.yml
# conda env export --no-builds -n rnaseq_tools > rnaseq_tools_env_full.yml

# To create environment.yml (cleaner version)
# conda env export --from-history -n isoseq_tools > isoseq_tools_env.yml
# conda env export --from-history -n rnaseq_tools > rnaseq_tools_env.yml
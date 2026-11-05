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
# Activate using: conda activate isoseq_tools 
# Deactivate using: conda deactivate


## SQANTI3 environment
#   This is part of the Isoseq analysis pipeline, but sqanti requires a separate conda env
wget https://github.com/ConesaLab/SQANTI3/releases/tag/v6.0/SQANTI3_v6.0.zip
mkdir sqanti3
unzip SQANTI3_v6.0.zip -d sqanti3
conda env create -f SQANTI3.conda_env.yml
# Activate using: conda activate sqanti3
# Deactivate using: conda deactivate


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
    subread \
    seqtk
# Activate using: conda activate rnaseq_tools 
# Deactivate using: conda deactivate


# To create environment.yml (with software versions)
# conda env export --no-builds -n isoseq_tools > isoseq_tools_env_full.yml
# conda env export --no-builds -n rnaseq_tools > rnaseq_tools_env_full.yml

# To create environment.yml (cleaner version)
# conda env export --from-history -n isoseq_tools > isoseq_tools_env.yml
# conda env export --from-history -n rnaseq_tools > rnaseq_tools_env.yml
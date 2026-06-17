#!/usr/bin/env bash
set -euo pipefail

## set up conda environment for Modern IsoSeq and RNA-Seq pipelines 

## IsoSeq environment 
#   Based on current IsoSeq CLI workflow (v4+)
#   No Cupcake, Python 3 only
conda create -n isoseq_tools -y \
  -c conda-forge \
  -c bioconda \
  python=3.10 \
  isoseq \
  lima \
  pbccs \
  pbmm2 \
  pbpigeon \
  pbbam \
  samtools \
  bamtools 


## RNASeq environment 
conda create -n rnaseq_tools -y \
  -c conda-forge \
  -c bioconda \
  python=3.10 \
  fastqc \
  star \
  trim-galore \
  subread \
  samtools 



# To create environment.yml 
# conda env export --no-builds -n isoseq_tools > isoseq_tools_env.yml
# conda env export --no-builds -n rnaseq_tools > rnaseq_tools_env.yml
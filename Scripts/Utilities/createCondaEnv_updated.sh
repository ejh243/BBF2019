#!/usr/bin/env bash
set -euo pipefail

## set up conda environment for Modern IsoSeq pipeline 
# Based on current IsoSeq CLI workflow (v4+)
# No Cupcake, Python 3 only

module load Miniconda3
source ./Config/config.txt

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


## set up conda environment for latest RNA-Seq pipeline 
# Latest version of FastQC and STAR already available on server 

module load FastQC
module load STAR

conda create -n rnaseq_tools -y \
  -c conda-forge \
  -c bioconda \
  trim-galore \
  subread


# To create environment.yml 
# conda env export --no-builds -n isoseq_tools > isoseq_tools_env.yml
# conda env export --no-builds -n rnaseq_tools > rnaseq_tools_env.yml
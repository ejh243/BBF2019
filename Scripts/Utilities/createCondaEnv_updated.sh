## set up conda environment Modern IsoSeq pipeline 
# Based on current IsoSeq CLI workflow (v4+)
# No Cupcake, Python 3 only

module load Miniconda3
source ./Config/config.txt

conda create -n isoseq_tools -y \
  -c conda-forge -c bioconda \
  python=3.10 \
  isoseq \
  lima \
  pbccs \
  pbmm2 \
  pigeon \
  samtools \
  bamtools \
  pysam

conda activate isoseq_tools


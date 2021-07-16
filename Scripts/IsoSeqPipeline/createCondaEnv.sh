## set up conda environment for isoseq3 pipeline
## although tutorial for install says python 3.7, lima only works with python 2.7

module load Miniconda2

conda create --name isoseq python=2.7
source activate isoseq

conda install -n isoseq biopython
conda install -n isoseq -c http://conda.anaconda.org/cgat bx-python

conda install -n isoseq -c bioconda isoseq3
conda install -n isoseq -c bioconda pbccs
conda install -n isoseq -c bioconda lima

## below are optional
conda install -n isoseq -c bioconda pbcoretools # for manipulating PacBio datasets
conda install -n isoseq -c bioconda bamtools    # for converting BAM to fasta
conda install -n isoseq -c bioconda pysam       # for making CSV reports

#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=48:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under.
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e IsoSeqAnalysis.err # error file
#PBS -o IsoSeqAnalysis.log # output file
#PBS -d /gpfs/ts0/projects/Research_Project-193495/BBF2019/ # set working directory to


module load Anaconda3/5.2.0
source activate isoseq
source ./config.txt

cd ${DATADIR}/BBF2019/

sh RunIsoSeq3_1Pipeline.sh
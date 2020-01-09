## post processing of iso-seq data

## create fq of high quality isoforms
gunzip Cluster/polished.hq.fasta.gz
fa2fq.py Cluster/polished.hq.fasta

## align with minimap2
minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 \
   ${REFGENOME} Cluster/polished.hq.fastq \
   > Aligned/hq_isoforms.fastq.sam \
  2> Aligned/hq_isoforms.fastq.sam.log
  
## sort sam file
samtools sort Aligned/hq_isoforms.fastq.sam > Aligned/hq_isoforms.sorted.sam
  
## convert sam to bam
samtools view -bS Aligned/hq_isoforms.fastq.sam > Aligned/hq_isoforms.fastq.bam
samtools sort Aligned/hq_isoforms.fastq.bam > Aligned/hq_isoforms.fastq.sorted.bam
samtools index Aligned/hq_isoforms.fastq.sorted.bam

## collapse redundant isoforms
collapse_isoforms_by_sam.py --input Cluster/polished.hq.fastq --fq \
   -s Aligned/hq_isoforms.sorted.sam --dun-merge-5-shorter -o collapsed_isoforms
## post processing of iso-seq data

## create fq of high quality isoforms
gunzip Cluster/polished.hq.fasta.gz
fa2fq.py Cluster/polished.hq.fasta
gzip Cluster/polished.hq.fasta
gzip Cluster/polished.hq.fastq

## align with minimap2
minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 \
   ${REFGENOME} Cluster/polished.hq.fasta.gz \
   > Aligned/hq_isoforms.fasta.sam \
  2> Aligned/hq_isoforms.fasta.sam.log
  
minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 \
   ${REFGENOME} Cluster/polished.hq.fastq.gz \
   > Aligned/hq_isoforms.fastq.sam \
  2> Aligned/hq_isoforms.fastq.sam.log
  
## sort sam file
samtools sort Aligned/hq_isoforms.fastq.sam > Aligned/hq_isoforms.sorted.sam
samtools sort Aligned/hq_isoforms.fasta.sam > Aligned/hq_isoforms.fasta.sorted.sam  

sort -k 3,3 -k 4,4n hq_isoforms.fastq.sam > hq_isoforms.fastq.sorted.sam
sort -k 3,3 -k 4,4n hq_isoforms.fasta.sam > hq_isoforms.fasta.sorted.sam

## convert sam to bam
samtools view -bS Aligned/hq_isoforms.fastq.sam > Aligned/hq_isoforms.fastq.bam
samtools sort Aligned/hq_isoforms.fastq.bam > Aligned/hq_isoforms.fastq.sorted.bam
samtools index Aligned/hq_isoforms.fastq.sorted.bam

## collapse redundant isoforms
#collapse_isoforms_by_sam.py --input Cluster/polished.hq.fastq --fq \
   -s Aligned/hq_isoforms.fastq.sorted.sam --dun-merge-5-shorter -o collapsed_isoforms

collapse_isoforms_by_sam.py --input Cluster/polished.hq.fasta \
   -s Aligned/hq_isoforms.fasta.sorted.sam --dun-merge-5-shorter -o Filtering/collapsed_isoforms.fasta
   
## get count information
fa2fq.py Filtering/collapsed_isoforms.fasta.collapsed.rep.fa
get_abundance_post_collapse.py Filtering/collapsed_isoforms.fasta.collapsed Cluster/polished.cluster_report.csv

## filter by number of full length reads for support
filter_by_count.py --min_count 2 --dun_use_group_count Filtering/collapsed_isoforms.fasta.collapsed

## filter away 5' degraded isoforms
filter_away_subset.py Filtering/collapsed_isoforms.fasta.collapsed


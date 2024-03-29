## prior to merging to get master transcriptome use sqanti3 to qc sample-level transcriptomes
## requires short read data aligned to individual transcriptomes
## input is collapsed high quality isoforms

p=$1
basename=${p%.subreads.bam}

echo "Running SQANTI3 on" ${basename}
mkdir -p ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3	

python ${SOFTWAREPATH}/SQANTI3/sqanti3_qc.py ${ALIGNEDDIR}/Collapsed/${basename}/out.collapsed.filtered.gff    \
   ${GENCODEGTF} ${REFGENOME} -o $basename -d ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3/    \
   --cage_peak ${SQANTIDATA}/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed    \
   --polyA_motif_list ${SQANTIDATA}/polyA_motifs/mouse_and_human.polyA_motif.txt    \
   -fl ${ALIGNEDDIR}/Collapsed/${basename}/out.collapsed.filtered.abundance.txt    \
   --isoAnnotLite --saturation \
   --expression $(ls -m ${GENECOUNTPATH}/RSEM/PersonalTranscriptome/${basename}/*.isoforms.results | tr -d '\n')    \
   -c $(ls -m ${ALIGNEDDIR}/SMRTcells/${basename}/*SJ.out.tab | tr -d '\n') \
   --cpus 4 --report pdf 
  
			 
## filter gtf based on sqanti classification
## don't filter on coverage at this point 
python ${SOFTWAREPATH}/SQANTI3/sqanti3_filter.py rules	   \
   --gtf ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3/${basename}_corrected.gtf   \
   --isoforms ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3/${basename}_corrected.fasta   \
   -d ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3/   \
   -o ${basename}   \
   ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3/${basename}_classification.txt


python ${SOFTWAREPATH}/SQANTI3/sqanti3_filter.py ML  \
   --gtf ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3/${basename}_corrected.gtf  \
   --isoforms ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3/${basename}_corrected.fasta   \
   -d ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3/   \
   -o ${basename}   \
   ${ALIGNEDDIR}/Collapsed/${basename}/SQANTI3/${basename}_classification.txt

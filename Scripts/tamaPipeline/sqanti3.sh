
mkdir -p ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3	

python ${SOFTWAREPATH}/SQANTI3/sqanti3_qc.py ${MASTERTRANSCRIPTOME}/TAMA/pfc_merge_smrt_all_nRead2_filtFrag.gtf    \
   ${GENCODEGTF} ${REFGENOME} -o TAMAmerge -d ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3    \
   --cage_peak ${SQANTIDATA}/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed    \
   --polyA_motif_list ${SQANTIDATA}/polyA_motifs/mouse_and_human.polyA_motif.txt   \
   --isoAnnotLite --saturation \
   --expression $(ls -m ${GENECOUNTPATH}/RSEM/TAMAmerge/${rnaID}*.isoforms.results | tr -d '\n')    \
   -c $(ls -m ${ALIGNEDDIR}/TAMAmerge/${rnaID}*SJ.out.tab | tr -d '\n') \
   --cpus 4 --report pdf 
  
## filter gtf based on sqanti classification
## don't filter on coverage at this point nb editted default parameters in 
python ${SOFTWAREPATH}/SQANTI3/sqanti3_filter.py rules	   \
   --gtf ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/TAMAmerge_corrected.gtf   \
   --isoforms ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/TAMAmerge_corrected.fasta   \
   -d ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/   \
   -o TAMAmerge_filter   \
   ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/TAMAmerge_classification.txt


## for interest apply ML approach
python ${SOFTWAREPATH}/SQANTI3/sqanti3_filter.py ML  \
   --gtf ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/TAMAmerge_corrected.gtf  \
   --isoforms ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/TAMAmerge_corrected.fasta \
   ${MASTERTRANSCRIPTOME}/TAMA/SQANTI3/TAMAmerge_classification.txt

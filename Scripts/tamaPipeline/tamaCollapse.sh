

p=$1
basename=${p%.subreads.bam}
echo "Collapsing " ${basename}

## if more than 250000 reads, split prior to running

NREADS=$(samtools view -c ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.sorted.bam)
if ${NREADS} > 250000:
then
   python ${SOFTWAREPATH}/tama/tama_go/split_files/tama_mapped_sam_splitter.py ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.sorted.bam 2 ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.sorted
   for i in {1..2}
   do
	python ${SOFTWAREPATH}/tama/tama_collapse.py -s ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.sorted_$i.sam -f ${REFGENOME} -p ${ALIGNEDDIR}/TAMA/Collapsed/${basename}_$i -x no_cap -rm low_mem
   done
	
else

    python ${SOFTWAREPATH}/tama/tama_collapse.py -s ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.sorted.bam -f ${REFGENOME} -p ${ALIGNEDDIR}/TAMA/Collapsed/${basename} -x no_cap -b BAM -rm low_mem
fi


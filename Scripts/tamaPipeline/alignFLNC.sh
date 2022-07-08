
p=$1
basename=${p%.subreads.bam}
echo "Aligning " ${basename}

## convert to fa file for alignment
bamtools convert -format fasta -in ${PROCESSEDDIR}/Refine/${basename}.flnc.bam > ${PROCESSEDDIR}/Refine/${basename}.flnc.fasta

## align
minimap2 -ax splice -t 30 -uf --secondary=no -C5 ${REFGENOME} ${PROCESSEDDIR}/Refine/${basename}.flnc.fasta > ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.sam


## convert to bam
samtools view -S -b ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.sam > ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.bam

## sort
samtools sort ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.bam -o ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.sorted.bam

rm ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.sam
rm ${ALIGNEDDIR}/TAMA/${basename}.all_isoforms.bam
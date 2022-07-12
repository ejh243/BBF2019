## SQANTI requires short read data aligned to the isoseq transcriptome
## this script generates gene/isoform counts matrix using RSEM which internally generates the STAR index for personalized transcriptome and aligns RNA seq data to it 
## had to split into separate STAR and RSEM steps as STAR default --limitOutSJcollapsed parameter was too low.
## recommended to run each sample in turn rather than as a single sample
## find trimmed files

all_f1=($(ls ${RNASeqDIR}Trimmed/*[rR]1*f*z))
all_f2=($(ls ${RNASeqDIR}Trimmed/*[rR]2*f*z))


GFF=${MASTERTRANSCRIPTOME}/TAMA/pfc_merge_smrt_all_nRead2_filtFrag.gtf

## use RSEM to generate counts matrix
mkdir -p ${RSEMREFDIR}/TAMAmerge/
rsem-prepare-reference --gtf ${GFF} ${REFGENOME} ${RSEMREFDIR}/TAMAmerge/TAMAmerge

# run STAR command with bespoke parameters
STAR  --runThreadN 8  --runMode genomeGenerate  --genomeDir ${RSEMREFDIR}/TAMAmerge/  --genomeFastaFiles ${REFGENOME}  --sjdbGTFfile ${GFF}  --sjdbOverhang 100  --outFileNamePrefix TAMAmerge --limitSjdbInsertNsj 1091177


mkdir -p ${ALIGNEDDIR}/TAMAmerge   
mkdir -p ${GENECOUNTPATH}/RSEM/TAMAmerge
for i in "${!all_f1[@]}"
do
    star_f1=${all_f1[$i]}
    star_f2=${all_f2[$i]}
    rnaID=$(basename ${star_f1})
    rnaID=${rnaID%%_*}
    
    echo "Aligning short read data for sample " ${rnaID}
    
	if [ ! -f ${GENECOUNTPATH}/RSEM/TAMAmerge/${rnaID}.isoforms.results ]
    then
     
        ## align with STAR 
        STAR --genomeDir ${RSEMREFDIR}/TAMAmerge/ --runThreadN 18 --readFilesIn ${star_f1},${star_f2} \
            --readFilesCommand zcat \
            --outFileNamePrefix ${ALIGNEDDIR}/TAMAmerge/${rnaID} \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outSAMattributes  NH HI AS NM MD \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999  \
            --outFilterMismatchNoverLmax 0.04  \
            --alignIntronMin 20  \
            --alignIntronMax 1000000  \
            --alignMatesGapMax 1000000  \
            --alignSJoverhangMin 8  \
            --alignSJDBoverhangMin 1  \
            --sjdbScore 1  \
            --genomeLoad NoSharedMemory  \
            --quantMode TranscriptomeSAM  \
            --outSAMheaderHD \@HD VN:1.4 SO:unsorted \
            --limitOutSJcollapsed 2000000 
        
        
        
        rsem-calculate-expression --num-threads 10 --alignments ${ALIGNEDDIR}/TAMAmerge/${rnaID}Aligned.toTranscriptome.out.bam ${RSEMREFDIR}/TAMAmerge ${GENECOUNTPATH}/RSEM/TAMAmerge/${rnaID}
        rm ${ALIGNEDDIR}/TAMAmerge/${rnaID}Aligned.out.bam
    fi
done


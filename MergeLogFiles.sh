## record size of raw BAM files

ls -lh RawData/*.subreads.bam > RawData/RawDataTable.txt

## merge log files from CCS step

ccsReport=($(ls CCS/*_report.txt))

join -t : ${ccsReport[0]} ${ccsReport[1]} > CCS/merged.ccs_report.txt

for f in ${ccsReport[@]:2}
do
	join -t : CCS/merged.ccs_report.txt $f > tmp.txt
	mv tmp.txt CCS/merged.ccs_report.txt
done

rm tmp.txt


## merge log files from lima step

limaReport=($(ls Lima/*.lima.summary))

join -t : ${limaReport[0]} ${limaReport[1]} > Lima/merged.lima.summary

for f in ${limaReport[@]:2}
do
	join -t : Lima/merged.lima.summary $f > tmp.txt
	mv tmp.txt Lima/merged.lima.summary
done

rm tmp.txt

## merge output of refine step


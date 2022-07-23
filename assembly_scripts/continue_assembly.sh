raw_dir=data/raw_reads
R1=$raw_dir/5-129B_TCACCAGGAC-TCCTTGTCTC_L001_R1_001.fastq.gz
R2=$raw_dir/5-129B_TCACCAGGAC-TCCTTGTCTC_L001_R2_001.fastq.gz
result_dir=results/assembly

coronaspades.py \
-t 1 \
-k 21,49,55 \
--restart-from k49 \
-o $result_dir


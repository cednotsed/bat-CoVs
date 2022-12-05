for rank in S G F O
do
python parse_kraken2_report.py \
	--i ../results/classification_out/kraken2_reports/ \
	--o ../results/classification_out/abundance_matrix/ \
	--rank $rank \
	--prefix abundance_matrix \
	--delim .

done

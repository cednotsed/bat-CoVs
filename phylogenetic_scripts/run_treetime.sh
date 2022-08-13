#prefix=all_deer_090722.n318.parsed_names.unambiguous_dates.v8_masked
#aln_file=data/alignments/${prefix}.aln
#tree_file=data/trees/${prefix}.contree
#date_file=data/metadata/all_deer_090722.n318.parsed.unambiguous_dates.dates_only.tsv
#res_dir=results/timetree/$prefix
#
#ls $aln_file
#ls $tree_file
#ls $date_file
#
#treetime \
#	--tree ${tree_file} \
#	--aln ${aln_file} \
#	--dates ${date_file} \
#	--clock-filter 3 \
#	--covariation \
#	--confidence \
#	--coalescent skyline \
#	--relax 1.0 0 \
#	--outdir ${res_dir}

#treetime clock \
#--tree ${tree_file} \
#--dates ${date_file} \
#--clock-filter 3 \
#--aln ${aln_file} \
#--relax 1.0 0 \
#--outdir ${res_dir}

aln_dir=i../data/alignments

for aln_file in $aln_dir/*pos.aln
do
	prefix=$(echo $aln_file| sed "s|$aln_dir/||g"| sed "s|.aln||g")
	date_file=data/metadata/${prefix}.v8_masked.dates_only.tsv
	tree_file=data/trees/${prefix}.treefile
	res_dir=results/timetree/$prefix
	echo $prefix
	ls $aln_file
	ls $date_file
	ls $tree_file
	echo $res_dir

	treetime \
		--tree ${tree_file} \
		--aln ${aln_file} \
		--dates ${date_file} \
		--clock-filter 3 \
		--covariation \
		--confidence \
		--coalescent skyline \
		--relax 1.0 0 \
		--keep-root \
		--outdir ${res_dir}

done

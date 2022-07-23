wkdir=../results/assembly/constructed_genomes
out_dir=$wkdir/paf_files
ref_dir=$wkdir/references
contig_dir=$wkdir/contig_pools

for contig_pool in $contig_dir/*
do
	ref=$ref_dir/$(echo $contig_pool|cut -d'_' -f5)
	sample_id=$(echo $contig_pool|cut -d'/' -f6 | cut -d'_' -f1)
	ls $ref
	echo $sample_id
	ls $out_dir
	minimap2 \
		-B 1 \
		-A 20 \
		$ref \
		$contig_pool \
		-o $out_dir/${sample_id}.tsv
done

prefix=coronaviridae.n2089_subset
prokka_out=../data/proteins/prokka_proteins/$prefix
out_dir=../results/pangenome/${prefix}_panaroo
input=$prokka_out/*.gff
ls $input
ls $out_dir

panaroo \
	-i $input \
	-o $out_dir \
	--clean-mode strict \
	-t 8 \
	-a core \
	--aligner mafft \
	--core_threshold 0.95 \
	-f 0.3 \
	--merge_paralogs \
	--refind_prop_match 0.5 \
	--search_radius 1000





prefix=coronaviridae.n2089_subset
prokka_out=../data/proteins/prokka_proteins/$prefix
out_dir=../results/pangenome/$prefix
input=$prokka_out/*.gff
ls $input
ls $out_dir

roary \
	-i 30 \
	-p 8 \
	-cd 90 \
	-s \
	-ap \
	-v \
	-e \
	-n \
	-f $out_dir \
	$input





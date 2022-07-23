db=./checkv-db-v1.2
db=../databases/checkv-db-v1.2/
assembly_dir=../results/assembly/
out_dir=../results/assembly/checkv_out

for i in $assembly_dir/all_novel*
do
	prefix=$(echo $i| sed "s|$assembly_dir/||g"| sed "s|_scaffolds.fasta||g")
    prefix=$(echo $i| sed "s|$assembly_dir/||g"| sed "s|_genomes.fna||g")
	contigs=$i
	sub_dir=$out_dir/$prefix
	echo $prefix
	echo $sub_dir
	#mkdir $sub_dir

	checkv end_to_end \
		-t 1 \
		-d $db \
		--remove_tmp \
		--restart \
		$contigs $sub_dir

done

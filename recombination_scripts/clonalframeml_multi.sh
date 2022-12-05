newick_dir=../data/trees/for_recombination_analysis
seq_dir=../data/alignments/for_recombination_analysis
out_dir=../results/recombination_out

for newick_file in $newick_dir/*.tree
do
	seq_file=$(echo $newick_file| sed "s|$newick_dir|$seq_dir|g"| sed "s|\\.tree|\\.aln|g")
	output_prefix=$(echo $newick_file| sed "s|$newick_dir|$out_dir|g"| sed "s|\\.tree||g")
	echo $newick_file
	echo $seq_file
	echo $out_prefix

	ClonalFrameML \
   		$newick_file $seq_file $output_prefix \
   		-emsim 100
done
#ClonalFrameML \
#	$newick_file $seq_file $output_prefix \
#	-emsim 100

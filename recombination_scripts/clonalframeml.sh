newick_file=../data/trees/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.contree
seq_file=../data/alignments/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.aln
output_prefix=../results/recombination_out/Sarbecovirus.trim_to_28261pos.n218
output_prefix=$(echo $seq_file|sed "s|data/alignments/for_recombination_analysis/|results/recombination_out/|g"| sed "s|.aln||g")

echo $output_prefix

ClonalFrameML \
	$newick_file $seq_file $output_prefix \
	-emsim 100

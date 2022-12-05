#newick_file=../data/trees/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.contree
#seq_file=../data/alignments/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.aln
#output_prefix=../results/recombination_out/gubbins_out/Sarbecovirus.trim_to_28261pos.n218

seq_file=../data/alignments/for_recombination_analysis/Sarbecovirus.n218.masked.aln
output_prefix=../results/recombination_out/gubbins_out/Sarbecovirus.n218.masked
run_gubbins.py \
	--prefix $output_prefix \
	--extensive-search \
	--threads 10 \
	--tree-builder iqtree-fast \
	-b 500 \
	$seq_file

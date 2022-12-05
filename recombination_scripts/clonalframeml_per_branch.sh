#newick_file=../data/trees/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.contree
#seq_file=../data/alignments/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.aln
#output_prefix=../results/recombination_out/clonalframeml_per_branch/Sarbecovirus.trim_to_28261pos.n218

newick_file=../data/trees/for_recombination_analysis/Sarbecovirus.n218.masked.tree.contree
seq_file=../data/alignments/for_recombination_analysis/Sarbecovirus.n218.masked.aln
output_prefix=../results/recombination_out/clonalframeml_per_branch/Sarbecovirus.n218.masked

ClonalFrameML \
    $newick_file $seq_file $output_prefix \
	-embranch true \
    -emsim 100 \
	-embranch_dispersion 0.1 \
#	-initial_values "0.076 0.00208 0.063"

newick_file=../data/trees/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.contree
seq_file=../data/alignments/for_recombination_analysis/Sarbecovirus.trim_to_28261pos.n218.aln
output_prefix=../results/recombination_out/Sarbecovirus.trim_to_28261pos.n218
output_prefix=../results/recombination_out/clonalframeml_per_branch/Sarbecovirus.trim_to_28261pos.n218

for emsim_file in ../results/recombination_out/*.emsim.txt
do
	output_prefix=$(echo $emsim_file|sed "s|.emsim.txt||g")
	Rscript cfml_results.R $output_prefix

done


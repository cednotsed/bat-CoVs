threads=10
#aln_path=../data/alignments/Sarbecovirus_subtree.n534.trim_to_28261pos.aln
#aln_path=../data/alignments/Pedacovirus_subtree.n106.trim_to_26674pos.aln
#aln_path=../data/alignments/Merbecovirus_subtree.n113.trim_to_26414pos.aln
aln_path=../data/alignments/for_recombination_analysis/Sarbecovirus.n218.masked.aln

output=$(echo $aln_path| sed "s|alignments|trees|g"| sed "s|.aln|.tree|g")
echo $output

iqtree2 \
	-nt $threads \
	-s $aln_path \
	-m GTR+G \
	--prefix $output \
	-B 1000 \
	-alrt 1000

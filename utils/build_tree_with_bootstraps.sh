threads=10
aln_path=../data/alignments/sarbecovirus.050422.trimmed.aln
output=$(echo $aln_path| sed "s|alignments|trees|g"| sed "s|.aln|.tree|g")
echo $output

iqtree2 \
	-nt $threads \
	-s $aln_path \
	-m GTR+G \
	--prefix $output \
	-B 1000 \
	-alrt 1000

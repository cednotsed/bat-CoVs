threads=10
reference=../data/genomes/references/NC_025217_sarbecovirus_outgroup.fna
genomes=../data/genomes/Sarbecovirus_subtree.n534.fna

reference=../data/genomes/references/MT663548_pedacovirus_outgroup.fna
genomes=../data/genomes/Pedacovirus_subtree.n106.fna

reference=../data/genomes/references/NC_009021_merbecovirus_outgroup.fna
genomes=../data/genomes/Merbecovirus_subtree.n113.fna

alignment_out=$(echo $genomes|sed 's/.fna/.aln/g'|sed "s|genomes|alignments|g")

echo $reference
echo $genomes
echo $alignment_out

augur align \
--nthreads ${threads} \
--sequences ${genomes} \
--output ${alignment_out} \
--reference-sequence ${reference} \
--method mafft


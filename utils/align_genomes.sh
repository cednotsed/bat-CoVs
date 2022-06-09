threads=10
reference=../data/genomes/references/bat_CoV_HKU9_NC_009021.1.fasta
genomes=../data/genomes/sarbecovirus.050422.fasta

reference=../data/genomes/references/WIV04_EPI_ISL_402124.fasta
genomes=../data/genomes/sarbecovirus_like.fasta

#reference=../data/genomes/references/MERS.ref.fasta
#genomes=../data/genomes/MERS_like.fasta

alignment_out=$(echo $genomes|sed 's/.fasta/.aln/g'|sed "s|genomes|alignments|g")

echo $reference
echo $genomes
echo $alignment_out

augur align \
--nthreads ${threads} \
--sequences ${genomes} \
--output ${alignment_out} \
--reference-sequence ${reference} \
--method mafft


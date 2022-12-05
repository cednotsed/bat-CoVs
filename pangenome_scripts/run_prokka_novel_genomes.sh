fasta_dir=../results/assembly/constructed_genomes
out_basedir=../results/prokka_out

# For 5-129B
fasta_dir=../results/assembly/constructed_genomes/5-129B_candidates

echo $out_basedir

for fasta in $fasta_dir/*.fna
do
	prefix=$(echo $fasta|sed "s|$fasta_dir/||g"|sed "s|.fna||g")
	echo $prefix
	echo $fasta
	prokka \
		--prefix $prefix \
		--kingdom Viruses \
		--evalue 1e-09 \
		--coverage 80 \
		--outdir $out_basedir \
		--force \
		$fasta
done


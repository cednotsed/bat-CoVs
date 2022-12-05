fasta_dir=../data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset
out_basedir=../results/prokka_out

echo $out_basedir

for fasta in $fasta_dir/NC_045512*.fna
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


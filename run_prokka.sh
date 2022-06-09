fasta_dir=data/genomes/raw_genomes/Alpha
fasta_dir=data/genomes/raw_genomes/Beta
out_basedir=data/proteins/prokka_proteins
protein_fasta=data/proteins/batch_entrez_proteins/human_infecting_CoVs.faa

for fasta in $fasta_dir/*.fasta
do
	echo $fasta
	prefix=$(echo $fasta|sed "s|$fasta_dir||g"|sed "s|/||g"|sed "s|\\.fasta||g")
	out_dir=$out_basedir/$prefix
	echo $out_dir
	prokka \
		--proteins $protein_fasta \
		--kingdom Viruses \
		--evalue 1e-09 \
		--coverage 80 \
		--metagenome \
		--outdir $out_dir \
		--force \
		$fasta
done

protein_fasta=data/proteins/batch_entrez_proteins/human_infecting_CoVs.faa

#fasta_dir=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/constructed_genomes
fasta_dir=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/data/genomes/human_infecting_CoVs

out_basedir=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/data/proteins/prokka_proteins/


for fasta in $fasta_dir/*.fna
do
	echo $fasta
	prefix=$(echo $fasta|sed "s|$fasta_dir||g"|sed "s|/||g"|sed "s|.fna||g"| sed "s|_assembly||g")
	out_dir=$out_basedir/$prefix
	echo $out_dir

	prokka \
		--proteins $protein_fasta \
		--kingdom Viruses \
		--evalue 1e-09 \
		--metagenome \
		--outdir $out_dir \
		--force \
		$fasta
done

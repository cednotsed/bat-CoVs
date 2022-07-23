genome_dir=./data/genomes

fasta=$genome_dir/alpha_merged.050422.fasta
#fasta=$genome_dir/beta_merged.050422.fasta

out_dir=./results/mash_out

out_path=$(echo $fasta|sed "s|\\.fasta||g"| sed "s|$genome_dir|$out_dir|g")

echo $out_path

mash sketch \
-i $fasta \
-o $out_path \
-k 12

mash dist ${out_path}.msh ${out_path}.msh -t > ${out_path}.tsv

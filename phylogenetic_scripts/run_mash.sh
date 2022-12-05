genome_dir=../data/genomes

fasta=$genome_dir/coronaviridae_n2118.080822.with_novel.fna

out_dir=../results/phylogenetic_out/mash_out

out_path=$(echo $fasta|sed "s|\\.fasta||g"| sed "s|$genome_dir|$out_dir|g")

echo $out_path

mash sketch \
-i $fasta \
-o $out_path \
-k 12

mash dist ${out_path}.msh ${out_path}.msh -t > ${out_path}.tsv

query_dir=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/constructed_genomes
ls $query_dir/*.fna > ./queries.txt
ref_dir=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/constructed_genomes/references
ls $ref_dir/*.fna > ./refs.txt

out=../results/phylogenetic_out/fastANI_out/average_nuc_ident.tsv

fastANI --ql queries.txt --rl refs.txt -o $out

ref=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/constructed_genomes/novel_genomes_to_send/all_novel_genomes.n9.fna
db=../databases/novel_genomes.n9.blastn_db/novel_genomes.n9.blastn_db

makeblastdb -input_type fasta \
    -in $ref \
    -out $db \
    -parse_seqids \
    -title "novel_genomes.n9.blastn_db" \
    -dbtype nucl

# Blastn
out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_primers_against_novel_genomes.tsv
query=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/data/genomes/primers.fna

echo ${query}
echo ${out}
echo ${db}

~/ncbi-blast-2.11.0+/bin/blastn -db ${db} \
	-task blastn-short \
   	-query ${query} \
   	-out ${out} \
	-qcov_hsp_perc 99 \
    -outfmt 6 \
    -num_threads 8


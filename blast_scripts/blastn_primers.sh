#ref=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/constructed_genomes/novel_genomes_to_send/all_novel_genomes.n9.fna
#db=../databases/novel_genomes.n9.blastn_db/novel_genomes.n9.blastn_db
#
#makeblastdb -input_type fasta \
#    -in $ref \
#    -out $db \
#    -parse_seqids \
#    -title "novel_genomes.n9.blastn_db" \
#    -dbtype nucl

#out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_primers_against_novel_genomes.tsv

# For blasting coronavirus database
ref=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/data/genomes/coronaviridae_n2118.080822.with_novel.fna
db=../databases/coronaviridae.n2118_subset.with_novel.blastn_db/coronaviridae.n2118_subset.with_novel.blastn_db
makeblastdb -input_type fasta \
    -in $ref \
    -out $db \
    -parse_seqids \
    -title "coronaviridae.n2118_subset.with_novel.blastn_db" \
    -dbtype nucl

out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_primers_against_all_genomes.tsv

# Blastn
query=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/data/genomes/all_primers.fna

echo ${query}
echo ${out}
echo ${db}

~/ncbi-blast-2.11.0+/bin/blastn -db ${db} \
	-task blastn-short \
   	-query ${query} \
   	-out ${out} \
	-qcov_hsp_perc 90 \
	-num_alignments 1000000000 \
	-evalue 1 \
    -outfmt 6 \
    -num_threads 8


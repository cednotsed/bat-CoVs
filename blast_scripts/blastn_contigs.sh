## Make db ##
#ref=../data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset.fna
#db=../databases/coronaviridae.n2089_subset.blastn_db/coronaviridae.n2089_subset.blastn_db
#makeblastdb -input_type fasta \
#	-in $ref \
#	-out $db \
#	-parse_seqids \
#	-title "coronaviridae.n2089_subset" \
#	-dbtype nucl

#ref=../data/genomes/UK_bats/all_UK_bat_sequences.fna
#db=../databases/UK_bats.blastn_db/UK_bats.blastn_db
#makeblastdb -input_type fasta \
#    -in $ref \
#    -out $db \
#    -parse_seqids \
#    -title "UK_bats.blastn_db" \
#    -dbtype nucl

#ref=../data/genomes/UK_bats/receptor_genes/bat_ACE2.fna
#db=../databases/bat_ACE2.blastn_db/bat_ACE2.blastn_db
#makeblastdb -input_type fasta \
#    -in $ref \
#    -out $db \
#    -parse_seqids \
#    -title "bat_ACE2.blastn_db" \
#    -dbtype nucl

#ref=../data/genomes/UK_bats/receptor_genes/bat_DPP4.fna
db=../databases/bat_DPP4.blastn_db/bat_DPP4.blastn_db
#makeblastdb -input_type fasta \
#    -in $ref \
#    -out $db \
#    -parse_seqids \
#    -title "bat_DPP4.blastn_db" \
#    -dbtype nucl

# Blastn
#out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_contigs.full_aln_only.tsv
#query=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/pooled_contigs.fna

#out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_all_novel_genomes.tsv
#query=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/all_novel_genomes.fna

#out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_to_UK_bats.tsv
#query=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/pooled_contigs.fna

out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_contigs_to_bat_DPP4.tsv
query=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/Plecotus_auritus_contigs.fna

#out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_contigs_to_bat_ACE2.tsv
#query=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/Rhinolophus_hipposideros_contigs.fna

echo ${query}
echo ${out}
echo ${db}

~/ncbi-blast-2.11.0+/bin/blastn -db ${db} \
	-task dc-megablast \
   	-query ${query} \
   	-out ${out} \
	-qcov_hsp_perc 80 \
    -outfmt 6 \
    -num_threads 8


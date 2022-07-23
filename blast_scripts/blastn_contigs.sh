## Make db ##
ref=../data/genomes/coronaviridae_taxid11118_complete_exclude_provirus_040722/coronaviridae.n2089_subset.fna
db=../databases/coronaviridae.n2089_subset.blastn_db/coronaviridae.n2089_subset.blastn_db
makeblastdb -input_type fasta \
	-in $ref \
	-out $db \
	-parse_seqids \
	-title "coronaviridae.n2089_subset" \
	-dbtype nucl

# Blastn
out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_contigs.full_aln_only.tsv
query=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/pooled_contigs.fna

out=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/blast_out/blast_all_novel_genomes.tsv
query=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/results/assembly/all_novel_genomes.fna

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


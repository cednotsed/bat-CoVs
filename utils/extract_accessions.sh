genome_file=../data/genomes/beta_merged.050422.fasta
accessions=../data/metadata/sarbecovirus.accessions_only.csv
output_file=../data/genomes/sarbecovirus.050422.fasta

seqtk subseq $genome_file $accessions > $output_file

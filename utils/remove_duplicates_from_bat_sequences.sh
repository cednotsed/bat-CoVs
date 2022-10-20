in_dir=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/data/genomes/UK_bats/UK_bats_raw
out_dir=/mnt/c/Users/Cedric/Desktop/git_repos/bat-CoVs/data/genomes/UK_bats/UK_bats_dedup

for i in $in_dir/*
do
	echo $i
	out_path=$(echo $i|sed "s|.fna|.dedup.fna|g"|sed "s|$in_dir|$out_dir|g")
	echo $out_path
	seqkit rmdup -s < $i > $out_path
done

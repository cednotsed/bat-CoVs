#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=120:0:0
#$ -wd /SAN/ballouxlab/uk_bats_meta/bat-CoVs/logs
#$ -S /bin/bash
#$ -j y
#$ -pe smp 24
#$ -R y
#$ -N kraken
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2

wkdir=/SAN/ballouxlab/uk_bats_meta/bat-CoVs
fastq_dir=/SAN/ballouxlab/uk_bats_meta_raw/raw_reads
out_dir=$wkdir/results/classification_out
tmp_dir=$wkdir/classification_scripts/kraken2_temp
n_threads=24
db=$wkdir/databases/k2_viral_20220607

#Classification
for fq1 in $fastq_dir/*R1*
do
	fq2=$(echo $fq1| sed "s|R1|R2|g")
    prefix=$(echo $fq1|sed "s|$fastq_dir/||g"| cut -d'_' -f1)
    echo $prefix
	out_path=$out_dir/$prefix.tsv

    kraken2 \
        --threads $n_threads\
        --db $db \
        --paired \
        --report $out_path \
        --output $tmp_dir/temp_output \
        $fq1 $fq2
done


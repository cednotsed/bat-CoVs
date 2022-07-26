#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -l h_rt=120:0:0
#$ -wd /SAN/ballouxlab/uk_bats_meta/bat-CoVs/logs
#$ -S /bin/bash
#$ -j y
#$ -pe smp 4
#$ -R y
#$ -N subset_reads
source ~/miniconda3/etc/profile.d/conda.sh
conda activate philodendron

in_dir=/SAN/ballouxlab/uk_bats_meta_raw/raw_reads
out_dir=/SAN/ballouxlab/uk_bats_meta_raw/subset_reads
n_reads=10000000
mkdir $out_dir

for fq1 in $in_dir/*R1*
do
    echo $fq1
    fq2=$(echo $fq1|sed "s|R1|R2|g")

    out1=$(echo $fq1|sed "s|$in_dir|$out_dir|g"|sed "s|001.|001.subset.|g")
    out2=$(echo $fq2|sed "s|$in_dir|$out_dir|g"|sed "s|001.|001.subset.|g")

    seqtk sample -s100 $fq1 $n_reads > $out1
    seqtk sample -s100 $fq2 $n_reads > $out2
done


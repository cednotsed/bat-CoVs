#$ -l tmem=2G
#$ -l h_vmem=2G
#$ -l h_rt=240:0:0
#$ -wd /SAN/ballouxlab/uk_bats_meta/bat-CoVs/logs
#$ -S /bin/bash
#$ -j y
#$ -pe smp 4
#$ -R y
#$ -N map_to_assembly
source ~/miniconda3/etc/profile.d/conda.sh
conda activate philodendron

bam_dir=/SAN/ballouxlab/uk_bats_meta/bat-CoVs/results/mapping_out/bam_files

for bam_path in $bam_dir/*.bam
do
    out_path=$(echo $bam_path|sed "s|bam_files|coverage_files|g"| sed "s|.bam|_coverage.tsv|g")
    echo $out_path
    samtools depth -a $bam_path > $out_path
done

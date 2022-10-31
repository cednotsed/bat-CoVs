#$ -l tmem=2G
#$ -l h_vmem=2G
#$ -l h_rt=120:0:0
#$ -wd /SAN/ballouxlab/uk_bats_meta/bat-CoVs/logs
#$ -S /bin/bash
#$ -j y
#$ -pe smp 8
#$ -R y
#$ -N Sample_5-129B
prefix=5-129B

source ~/miniconda3/etc/profile.d/conda.sh
conda activate philodendron

wkdir=/SAN/ballouxlab/uk_bats_meta/bat-CoVs
fastq_dir=/SAN/ballouxlab/uk_bats_meta_raw/subset_reads/kraken2_positive
bam_dir=$wkdir/results/classification_out/bam_files
flag_dir=$wkdir/results/classification_out/flagstat_files
cov_dir=$wkdir/results/classification_out/coverage_files
ref=$wkdir/data/genomes/references/MW719567.fna
n_threads=8

fq1=$fastq_dir/$prefix*R1*
fq2=$fastq_dir/$prefix*R2*
echo $prefix
bam_path=$bam_dir/$prefix.bam
flag_path=$flag_dir/$prefix.flagstat.txt
cov_path=$cov_dir/$prefix.coverage.txt

# Map
bwa mem -t $n_threads $ref $fq1 $fq2 | samtools view -bS -@ $n_threads - | samtools sort -@ $n_threads -n -o $bam_path -
samtools index $bam_path

# Flagstat
samtools flagstat $bam_path > $flag_path

# Coverage
samtools depth -@ $n_threads $bam_path > $cov_path


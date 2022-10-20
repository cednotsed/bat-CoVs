#$ -l tmem=3G
#$ -l h_vmem=3G
#$ -l h_rt=120:0:0
#$ -wd /SAN/ballouxlab/uk_bats_meta/bat-CoVs/logs
#$ -S /bin/bash
#$ -j y
#$ -pe smp 7
#$ -R y
#$ -N qc_4-126A
prefix=4-126A
source ~/miniconda3/etc/profile.d/conda.sh
conda activate assembly
# Start time
start=`date +%s`

raw_dir=/SAN/ballouxlab/uk_bats_meta_raw/raw_reads
echo $prefix

R1=$(ls $raw_dir/${prefix}*R1*.gz)
R2=$(ls $raw_dir/${prefix}*R2*.gz)
out_dir=/SAN/ballouxlab/uk_bats_meta_raw/novel_QCed
out1=$(echo $R1|sed "s|$raw_dir|$out_dir|g"|sed "s|.fastq.gz|.QC.fastq.gz|g")
out2=$(echo $R2|sed "s|$raw_dir|$out_dir|g"|sed "s|.fastq.gz|.QC.fastq.gz|g")

adapter_file=/SAN/ballouxlab/uk_bats_meta/bat-CoVs/assembly_scripts/adapters.fa

bbduk.sh in1=$R1 in2=$R2 out1=$out1 out2=$out2 \
    ref=$adapter_file ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    qtrim=r trimq=10 \
    maq=10

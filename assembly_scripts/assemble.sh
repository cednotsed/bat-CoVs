#$ -l tmem=20G
#$ -l h_vmem=20G
#$ -l h_rt=120:0:0
#$ -wd /SAN/ballouxlab/uk_bats_meta/bat-CoVs/logs
#$ -S /bin/bash
#$ -j y
#$ -pe smp 8
#$ -R y
#$ -N n_1-GH087_QCed_n8_160g
prefix=1-GH087
source ~/miniconda3/etc/profile.d/conda.sh
conda activate assembly
# Start time
start=`date +%s`

wkdir=/SAN/ballouxlab/uk_bats_meta/bat-CoVs
#raw_dir=/SAN/ballouxlab/uk_bats_meta_raw/raw_reads/pcr_positive
raw_dir=/SAN/ballouxlab/uk_bats_meta_raw/novel_QCed
echo $prefix
R1=$raw_dir/${prefix}*R1*.gz
R2=$raw_dir/${prefix}*R2*.gz
#result_dir=$wkdir/results/coronaspades_n8_m250/complete_batch2/$prefix
result_dir=$wkdir/results/coronaspades_n8_m250/novel_QCed/$prefix

coronaspades.py \
--threads 8 \
--memory 160 \
-1 $R1 \
-2 $R2 \
-o $result_dir

# End time
end=`date +%s`
runtime=$(($(($end - $start)) / 3600))
echo $start
echo $end
echo $runtime

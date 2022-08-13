#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=240:0:0
#$ -wd /SAN/ballouxlab/uk_bats_meta/bat-CoVs/logs
#$ -S /bin/bash
#$ -j y
#$ -pe smp 8
#$ -R y
#$ -N map_to_assembly
source ~/miniconda3/etc/profile.d/conda.sh
conda activate philodendron

ref_dir=/SAN/ballouxlab/uk_bats_meta/bat-CoVs/results/assembly/constructed_genomes
fastq_dir=/SAN/ballouxlab/uk_bats_meta_raw/subset_reads
out_dir=/SAN/ballouxlab/uk_bats_meta/bat-CoVs/results/assembly/bam_files

for ref in $ref_dir/*.fna
do
    prefix=$(echo $ref|sed "s|$ref_dir/||g"|cut -d'_' -f1)
    fq1=$fastq_dir/${prefix}_*R1*
    fq2=$fastq_dir/${prefix}_*R2*
    ls $fq1
    out_path=$out_dir/$prefix.bam
    echo $out_path

    # Index
    bowtie2-build $ref $ref
    bowtie2 \
        -x $ref \
        -1 $fq1 \
        -2 $fq2 \
        | samtools view -bS \
        | samtools sort - -o $out_path

    samtools index $out_path
done

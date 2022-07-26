#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=120:0:0
#$ -wd /SAN/ballouxlab/uk_bats_meta/bat-CoVs/logs
#$ -S /bin/bash
#$ -j y
#$ -pe smp 8
#$ -R y
#$ -N make_krakenuniq_db
source ~/miniconda3/etc/profile.d/conda.sh
conda activate krakenuniq

#genome_dir=../databases/coronaviridae.n2089_subset
#DBNAME=../databases/coronaviridae.n2089_subset_krakenuniq_db

#for file in $genome_dir/*.fna
#do
#        kraken-build --add-to-library $file --db $DBNAME
#done

DBNAME=/SAN/ballouxlab/uk_bats_meta/bat-CoVs/databases/refseq_viral
mkdir $DBNAME

#krakenuniq-download --db $DBNAME taxonomy
#krakenuniq-download --db $DBNAME --dust refseq/viral/Any viral-neighbors
#kraken-build --download-taxonomy --db $DBNAME
kraken-build --download-library viral --db $DBNAME
kraken-build --build --db $DBNAME --threads 8

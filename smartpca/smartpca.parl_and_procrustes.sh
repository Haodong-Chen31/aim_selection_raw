#!/bin/bash

geno_dir=$1
code_dir=$2
pop=$3
parallel_threads=$4

###################### smartpca
#pop=poplist_1kg.txt
#geno_dir=/mnt/hdd1/HaodongChen/aim_260114_1kg_top1to2000_rmoutlier/data/tagSNP_union/
#code_dir=$(pwd)

#计算aim集的主成分
for n in $(seq 1 1 2000);do
	echo "bash smartpca.sh aim${n} ${pop} ${geno_dir} ${code_dir}"
done > smartpca.parl
cat smartpca.parl | parallel -j ${parallel_threads}
wait

mkdir -p pca_results
mv smartpca_aim* PCs_aim* plot_aim* pca_results/


###################### procrustes
# 在文件中手动修改运行范围 line 18 19 24
# 传参为原始pca的evec文件（建议绝对路径）
Rscript run_procrustes_analysis.R "$(pwd)/smartpca_modern_EA_qc_kin3.evec"

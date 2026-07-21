#!/bin/bash

# 原始全基因组数据 质控后的
origin_fn="/mnt/hdd1/HaodongChen/aim_260301_arraysample/language_tag_0to2k/data/modern_EA_qc_kin3"
pop="poplist_HO.txt"
parallel_threads=32 # plink和smartpca线程数

cat > get_aim_dataset.sh << 'EOF'
#!/bin/bash

origin_fn=$1
n=$2

tagSNP_snpid=/mnt/hdd1/HaodongChen/aim_260301_arraysample/language_tag_0to2k/data/tagSNP_union/tagSNP_union.top${n}.snpid

# 从质控后的全基因组数据中提aim子集
plink --bfile ${origin_fn} --extract ${tagSNP_snpid} --make-bed --out aim${n}
bash convert_plink2eigen.sh aim${n}
EOF

mkdir -p data
# 获得不同aim集的数据
for n in $(seq 1 1 2000);do
	echo "bash get_aim_dataset.sh ${origin_fn} ${n}"
done > get_aim_dataset.parl
cat get_aim_dataset.parl | parallel -j ${parallel_threads}
mv aim* data/
rm par*

# 跑smartpca和普鲁克分析
# 根据aim top集修改line 14的seq
# 以及run_procrustes_analysis.R的line18 19 24 和 line32aim pca的文件名前缀
bash smartpca.parl_and_procrustes.sh data . ${pop} ${parallel_threads}

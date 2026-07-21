#!/bin/bash

n=$1
#top SNP画曼哈顿 average和CV比较
head -1 combined_fst.txt > combined_fst_top${n}.txt
awk 'NR==FNR {snp_ids[$1]; next} $3 in snp_ids' ./tagSNP_union/tagSNP_union.top${n}.snpid combined_fst.txt >> combined_fst_top${n}.txt
# mean fst  CV
Rscript ../code/SNPselect_method1.r combined_fst_top${n}.txt 0.05 0.1 > ../log/SNPselect_method1_top${n}.log
#染色体上SNP统计
echo -e "chrom,total" > ./tagSNP_union/tagSNP_union.top${n}_number_snp.csv
cat ./tagSNP_union/tagSNP_union.top${n}.bim | awk '{print $1}' | sort -V | uniq -c | awk '{print $2","$1}' >> ./tagSNP_union/tagSNP_union.top${n}_number_snp.csv
#本体画图
#Rscript ../code/number_snp.R tagSNP_union.top${n}_number_snp.csv

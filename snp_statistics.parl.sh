#!/bin/bash

parallel_thread=$1

cd $(dirname $(pwd))/data
for n in 1 $(seq 1 1 2000);do
	echo "bash ../code/snp_statistics.sh ${n}"
	#top SNP画曼哈顿 average和CV比较
	#head -1 combined_fst.txt > combined_fst_top${n}.txt
	#awk 'NR==FNR {snp_ids[$1]; next} $3 in snp_ids' ./tagSNP_union/tagSNP_union.top${n}.snpid combined_fst.txt >> combined_fst_top${n}.txt
	#Rscript ../code/SNPselect_method1.r combined_fst_top${n}.txt 0.1 0.5 > ../log/SNPselect_method1_top${n}.log
	#染色体上SNP统计
	#echo -e "chrom,total" > ./tagSNP_union/tagSNP_union.top${n}_number_snp.csv
	#cat ./tagSNP_union/tagSNP_union.top${n}.bim | awk '{print $1}' | sort -V | uniq -c | awk '{print $2","$1}' >> ./tagSNP_union/tagSNP_union.top${n}_number_snp.csv
	#本体画图
	#Rscript ../code/number_snp.R tagSNP_union.top${n}_number_snp.csv
done > snp_statistics.parl
cat snp_statistics.parl | parallel -j ${parallel_thread}

mkdir combined_fst
mv combined_fst_* combined_fst/

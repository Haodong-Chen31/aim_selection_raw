#!/bin/sh

for n in 1 $(seq 5 5 500);do
	echo "sh extract_pca.sh ../data/tagSNP_union/tagSNP_union.top${n}_qc poplist.txt > ../log/pca_tagSNP_union.top${n}_qc.log"
done > pca.parl
cat pca.parl | parallel -j 16

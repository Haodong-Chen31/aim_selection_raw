#!/bin/bash

bfile_fn=$1
parallel_thread=$2

pop=$(dirname $(pwd))/data/pop
cp convert_plink2eigen.sh $(dirname $(pwd))/data
for n in $(seq 1 1 2000);do
        echo "bash step2_test4.sh ${pop} ${n} ${bfile_fn} > ../log/step2_test4_top${n}.log"
done > step2_test4.parl
cat step2_test4.parl | parallel -j ${parallel_thread}

mkdir -p $(dirname $(pwd))/data/tagSNP_union
cd $(dirname $(pwd))/data
#mv tagSNP_union.* tagSNP_union/
find . -name "tagSNP_union*" -type f -print0 | xargs -0 mv -t tagSNP_union/
#zip -r tagSNP_union_1kg_rmoutliers_top1_to_top10000.zip tagSNP_union/

#统计位点数
cat ../log/step2_test4_top*.log | grep tagSNP数量 | sort -V > ../aims_number.txt
sed -i 's/的/ 的/g' ../aims_number.txt
# cat step2_test4.parl2 | parallel -j 5 --joblog ../log/parallel2_step2_test4.log

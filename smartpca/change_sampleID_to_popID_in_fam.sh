#!/bin/sh

pop_sample=$1
fam=$2

awk 'BEGIN {OFS="\t"} 
     NR==FNR {map[$2] = $1; next}  # 读取pop_sample_all797.txt，建立映射：样本ID -> 群体
     $2 in map {$1 = map[$2]; print}  # 如果fam第二列在映射中，更新第一列
     !($2 in map) {print}  # 如果没有匹配，保留原样（或可以设置为默认值）' \
    ${pop_sample} \
    ${fam}

#pop_sample_all797.txt 第一列为群体名，第二列为样本名

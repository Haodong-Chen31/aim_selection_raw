#!/bin/bash

origin_fn=$1
n=$2

tagSNP_snpid=/mnt/hdd1/HaodongChen/aim_260114_1kg_top1to2000_rmoutlier/data/tagSNP_union/tagSNP_union.top${n}.snpid

# 从质控后的全基因组数据中提aim子集
plink --bfile ${origin_fn} --extract ${tagSNP_snpid} --make-bed --out aim${n}
bash convert_plink2eigen.sh aim${n}

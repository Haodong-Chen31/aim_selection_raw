#!/bin/sh

BASE_DIR=$(dirname $(pwd))
PCA_DIR=${BASE_DIR}/smartpca
DATA_DIR=${BASE_DIR}/data

mkdir ${PCA_DIR}
cp pca.parl.sh extract_pca.sh convert_eigen2plink.sh ${PCA_DIR}

FN="merged_1kg_sgdp.alldone"
cd ${DATA_DIR}
plink --bfile ${FN} --geno 0.1 --mind 0.1 --hwe 1e-3 --make-bed --out ${FN}_qc
sh convert_plink2eigen.sh ${FN}_qc

# pca 在本地可视化
cp ${FN}_qc.{geno,snp,ind} ${PCA_DIR}
cd ${PCA_DIR}
cat ${FN}_qc.ind | awk '{print $3}' | sort | uniq > poplist.txt
sh extract_pca.sh ${FN}_qc poplist.txt
sh pca.parl.sh

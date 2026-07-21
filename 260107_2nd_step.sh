#!/usr/bin/env bash

set -e # 如果任何命令失败，脚本将立即退出

#export TMPDIR="/mnt/hdd1/HaodongChen/tmp"
#expiort TMP="$TMPDIR"
#export TEMP="$TMPDIR"

# --- 配置 ---
FN="merged_1kg_EAS_autosomes"
THREADS=8
MEM="16G" # 根据你的机器内存调整
haploview_dir=/home/HaodongChen/software/

# --- 路径设置 ---
# 假设此脚本位于 'code' 目录
BASE_DIR=$(dirname $(pwd)) # 获取项目根目录
DATA_DIR="${BASE_DIR}/data"
LOG_DIR="${BASE_DIR}/log"
CODE_DIR="${BASE_DIR}/code"
PCA_DIR="${BASE_DIR}/smartpca"

# nb mcc效能验证
echo "==== 开始效能验证... ===="
mkdir -p ${BASE_DIR}/nb_mcc ; cd ${BASE_DIR}/nb_mcc
cp -r ${DATA_DIR}/tagSNP_union/ .

Rscript nb_mcc_functions.R > ${LOG_DIR}/nb_mcc.log 2>&1

# smartpca
echo "==== 开始跑TOP SNP的smartpca... ===="
mkdir -p ${PCA_DIR} ; cd ${PCA_DIR}
# 根据top数确定并行数
bash smartpca.parl.sh 128

echo "===== PIPELINE COMPLETED SUCCESSFULLY ====="


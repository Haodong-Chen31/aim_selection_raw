#!/usr/bin/env bash

set -e # 如果任何命令失败，脚本将立即退出

#export TMPDIR="/mnt/hdd1/HaodongChen/tmp"
#expiort TMP="$TMPDIR"
#export TEMP="$TMPDIR"

# --- 配置 ---
FN="modern_EA_qc_kin3"
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

# --- 主流程 ---
echo "===== STARTING ANALYSIS PIPELINE ====="

# 1. 准备样本和人群文件
echo "--> Step 1: Preparing sample and population files..."
mkdir -p "${LOG_DIR}"
cd "${DATA_DIR}"
if [ ! -f "all_EA_samples.txt" ] || [ ! -f "pop" ]; then
    #sh sample_ind.sh # 如果需要，运行此脚本
    cat ind_*.txt > all_EA_samples.txt
    find . -name 'ind_*.txt' -printf '%f\n' | sed -e 's/^ind_//' -e 's/\.txt$//' | sort > pop
    awk '{split(FILENAME, f, "[_.]"); print $1"\t"f[2]}' ind_*.txt > ind2pop.txt
else
    echo "Sample/pop files already exist, skipping generation."
fi

# 2. 高效处理VCF: 筛选、排序、重命名、索引
#echo "--> Step 2: High-performance VCF processing..."
#if [ ! -f "${FN}.vcf.gz" ]; then
#    bcftools view -S all_EA_samples.txt merged_autosome_qc.vcf.gz -Ou | \
#        bcftools sort -m ${MEM} -Ou | \
#        bcftools annotate --set-id '%CHROM\_%POS' -Oz -o ${FN}.vcf.gz --threads ${THREADS}
#    
#    bcftools index --threads ${THREADS} ${FN}.vcf.gz
#else
#    echo "Final VCF ${FN}.vcf.gz already exists, skipping generation."
#fi

# 3. 转换为PLINK格式并更新FAM和表型
# echo "--> Step 3: Converting to PLINK and updating files..."
# if [ ! -f "${FN}.bed" ]; then
#     plink --vcf ${FN}.alldone.rmoutliers.vcf.gz --make-bed --const-fid --out ${FN}
# 
#     # 更新Family ID为Population ID
#     awk 'BEGIN{OFS="\t"} NR==FNR {pop[$1]=$2; next} {$1=(($2 in pop) ? pop[$2] : "UNKNOWN"); print}' ind2pop.txt ${FN}.fam > tmp.fam
#     mv tmp.fam ${FN}.fam
# 
#     # 更新表型为1 (unaffected), 直接在二进制文件上操作
#     #plink --bfile ${FN} --make-bed --pheno-name 1 --out ${FN}_tmp
#     #mv ${FN}_tmp.bed ${FN}.bed
#     #mv ${FN}_tmp.bim ${FN}.bim
#     #mv ${FN}_tmp.fam ${FN}.fam
# else
#     echo "BPLINK fileset ${FN} already exists, skipping generation."
# fi

if [ ! -f "${FN}.ped" ]; then
    plink --bfile ${FN} --recode --out ${FN}

    # 更新Family ID为Population ID

    # 更新表型为1 (unaffected)
    cat ${FN}.ped | awk '$6=1{print}' > tmp.ped
    rm ${FN}.ped ; mv tmp.ped ${FN}.ped
else
    echo "PLINK fileset ${FN} already exists, skipping generation."
fi

# 4. 执行下游并行分析
echo "--> Step 4: Running downstream parallel analyses..."
cd "${CODE_DIR}"

#fst r2
echo "==== 开始计算fst... ===="
bash step1_aim_tagSNP_select.parl.sh ${FN} 500 0.05 0.1 ../data/pop 32 ${haploview_dir}

# 根据top数确定并行数
echo "==== 开始筛选位点... ===="
bash step2_test4.parl.sh ${FN} 128 | tee ../log/step2.log

echo "==== 开始合并和统计结果... ===="
bash combined_fst.sh ../data/pop ${FN} > ../log/combined_fst.log

# 根据top数确定并行数
bash snp_statistics.parl.sh 128

# nb mcc效能验证
echo "==== 开始效能验证... ===="
mkdir -p ${BASE_DIR}/nb_mcc ; cd ${BASE_DIR}/nb_mcc
cp -r ${DATA_DIR}/tagSNP_union/ .
cp ${CODE_DIR}/nb_mcc/* .
bash all_steps.sh

# 记得改参数
#Rscript nb_mcc_functions.R > ${LOG_DIR}/nb_mcc.log 2>&1

# smartpca
#echo "==== 开始跑TOP SNP的smartpca... ===="
mkdir -p ${PCA_DIR} ; cd ${PCA_DIR}
cp ${CODE_DIR}/smartpca/* .
# 根据top数确定并行数
#bash smartpca.parl.sh 128

#先跑一个原始的
bash smartpca.sh modern_EA_qc_kin3 poplist_HO.txt /mnt/hdd1/HaodongChen/aim_260301_arraysample/language_tag_0to2k/data .

#aim的pca以及普鲁克
bash smartpcaandprocrustes_using_aim.sh

echo "===== PIPELINE COMPLETED SUCCESSFULLY ====="


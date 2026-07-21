#!/bin/bash

#step2_test4.sh
#取每个种群对top number的SNP，并取并集
mapfile -t pops < $1
top_num=$2
bfile_fn=$3
tagSNP_pop1_pop2=""
missing_files=""
all_files_exist=true

cd $(dirname $(pwd))/data
for ((i=0; i<${#pops[@]}; i++)); do
    for ((j=i+1; j<${#pops[@]}; j++)); do
        # 定义当前种群对的文件名
        cat ./test_${pops[i]}_${pops[j]}/tagSNP_${pops[i]}_${pops[j]}.txt | sort -rk4,4 | head -${top_num} | awk '{print $1}' > ./test_${pops[i]}_${pops[j]}/tagSNP_${pops[i]}_${pops[j]}.top${top_num}.snpid
        top_tag_file="./test_${pops[i]}_${pops[j]}/tagSNP_${pops[i]}_${pops[j]}.top${top_num}.snpid"

        # 检查是否每组都生成了对应的文件
        if [[ -f ${top_tag_file} ]]; then
            # 将文件名添加到数组中
            tagSNP_pop1_pop2="${tagSNP_pop1_pop2} ${top_tag_file}"
            echo "找到文件: ${top_tag_file}"
        else
            echo "错误: 文件不存在 - ${top_tag_file}"
            all_files_exist=false
            missing_files="${missing_files} ${top_tag_file}"
        fi
    done
done

# 如果有任何文件缺失，则退出不执行后续代码
if [[ ${all_files_exist} != "true" ]]; then
    echo "错误: 以下文件缺失: ${missing_files}"
    echo "终止执行，不进行SNP交集和并集计算"
    exit 1
fi

# 所有文件都存在，执行SNP的交集和并集计算
echo "所有文件都存在，开始计算SNP的交集和并集..."
echo "找到 $(echo ${tagSNP_pop1_pop2} | wc -w) 个tagSNP.top${top_num}文件"

#并集
sort ${tagSNP_pop1_pop2} | uniq > tagSNP_union.top${top_num}.snpid
echo "top${top_num}的tagSNP数量（并集）： $(wc -l < tagSNP_union.top${top_num}.snpid)"
#获得snp位置文件
cat tagSNP_union.top${top_num}.snpid | sed 's/_/\t/g' | sort -V > tagSNP_union.top${top_num}.snpsite
# snpsite转BED
#awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1}' tagSNP_union.top${top_num}.snpsite > tagSNP_union.top${top_num}.BED
# liftOver转基因组坐标
#~/software/liftOver tagSNP_union.top${top_num}.BED /mnt/hdd7/HaodongChen/aim_260105/chain_files_human/GRCh38_to_GRCh37.chain.gz tagSNP_union.top${top_num}.liftover2grch37.BED tagSNP_union.top${top_num}.unmap.BED

#
plink --bfile ${bfile_fn} --extract tagSNP_union.top${top_num}.snpid --make-bed --out tagSNP_union.top${top_num}
# 转为dosage格式
#plink --bfile tagSNP_union.top${top_num} --recode A --out tagSNP_union.top${top_num}.dosage
#sed -i '1{s/_[ACGT]//g}' tagSNP_union.top${top_num}.dosage.raw

plink --bfile tagSNP_union.top${top_num} --geno 0.1 --mind 0.1 --hwe 1e-3 --make-bed --out tagSNP_union.top${top_num}_qc
bash convert_plink2eigen.sh tagSNP_union.top${top_num}_qc
rm parconvert_plink2eigen_tagSNP_union.top${top_num}_qc

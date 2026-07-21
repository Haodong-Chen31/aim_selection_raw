#!/usr/bin/env bash

# 检查参数数量
if [ $# -ne 7 ]; then
    echo "错误: 参数数量不正确!"
    echo "使用方法: $0 <prefix of input file> <ld_window_kb> <fst_threshold> <r2_threshold> <input population file> <parallel_threads> <haploview_dir>"
    echo "参数说明:"
    echo "  <prefix of input file>: PLINK二进制文件前缀 (.bed/.bim/.fam)"
    echo "  <ld_window_kb>: LD计算窗口大小 (单位: kb)"
    echo "  <fst_threshold>: FST统计量阈值，用于筛选高分化位点"
    echo "  <r2_threshold>: 连锁不平衡R²阈值，用于筛选独立位点"
    echo "  <input population file>: 包含种群列表的文本文件(!!含路径!!)，每行一个种群名称"
    echo "  <parallel_threads>: 并行处理的线程数"
    echo "  <haploview_dir>：haploview软件目录"
    exit 1
fi

fn_input=$1
ld_window_kb=$2
fst_threshold=$3
r2_threshold=$4
mapfile -t pops < $5
thread_parallel=$6
haploview_dir=$7
dirname=$(dirname $(pwd))

echo "开始标签SNP筛选流程..."
echo "输入文件: ${fn_input}"
echo "LD窗口: ${ld_window_kb}kb"
echo "FST阈值: ${fst_threshold}"
echo "R²阈值: ${r2_threshold}"
echo "种群数量: ${#pops[@]}"
echo "并行线程: ${thread_parallel}"

# mkdir -p ${dirname}/log
# #思路：计算pairwise population的fst；plink计算两两位点的r2连锁不平衡系数（fst>0.01,r2<0.1）
# #计算500kb物理距离内的SNP对，只输出r² ≥ 0.1的结果
# cd ${dirname}/data
# #改SNP的ID
# #zcat EA_samples_autosomes.vcf.gz | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$3 = $1 "_" $2; print}' | bgzip > EA_samples_autosomes_snp_renamed.vcf.gz
# #rm EA_samples_autosomes.vcf.gz ; mv EA_samples_autosomes_snp_renamed.vcf.gz EA_samples_autosomes.vcf.gz
# bcftools index EA_samples_autosomes.vcf.gz
# plink --vcf EA_samples_autosomes.vcf.gz --make-bed --out EA_samples_autosomes
# # sort
# #bcftools sort EA_samples_autosomes.vcf.gz -Oz -o EA_samples_autosomes.sorted.vcf.gz
# # 修改fam
# > ind2pop.txt
# for f in $(ls ind_*.txt);do pop=$(echo ${f} | sed 's/^ind_//g;s/.txt$//g') ; awk -v pop=${pop} '{print $1"\t"pop}' ${f} >> ind2pop.txt;done
# awk 'BEGIN{OFS="\t"}
#      NR==FNR {pop[$1]=$2; next}
#      {
#        if ($2 in pop) $1=pop[$2];
#        else $1="UNKNOWN";
#        print
#      }' ind2pop.txt EA_samples_autosomes.fam > tmp.fam
# rm EA_samples_autosomes.fam;mv tmp.fam EA_samples_autosomes.fam
# plink --bfile EA_samples_autosomes --recode --out EA_samples_autosomes
# cat EA_samples_autosomes.ped | awk '$6=1{print}' > tmp.ped
# rm EA_samples_autosomes.ped ; mv tmp.ped EA_samples_autosomes.ped

cd ${dirname}/data
plink --bfile ${fn_input} --r2 --ld-window-kb ${ld_window_kb} --ld-window-r2 ${r2_threshold} --out ${fn_input}_r2_morethan${r2_threshold}

#计算pairwise population的fst
#经过测试，vcf文件中有其他人群对计算成对人群的fst没有影响
#data/下存储了所有人群，以及每个人群的样本ID
for ((i=0; i<${#pops[@]}; i++)); do
    for ((j=i+1; j<${#pops[@]}; j++)); do
        pop1_pop2=${pops[i]}-${pops[j]}
        echo "bash ${dirname}/code/aim_tagSNP_select.sh ${pops[i]} ${pops[j]} ${fn_input} ${fst_threshold} ${r2_threshold} ${haploview_dir} > ${dirname}/log/aim_tagSNP_select_${pops[i]}_${pops[j]}.log"
    done
done > aim_tagSNP_select.parl

cat aim_tagSNP_select.parl | parallel -j ${thread_parallel} --progress --joblog aim_tagSNP_select.parl.log

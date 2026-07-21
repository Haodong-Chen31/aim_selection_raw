#!/usr/bin/env bash

pop1=$1
pop2=$2
fn=$3
fst_threshold=$4
r2_threshold=$5
haploview_dir=$6

cd $(dirname $(pwd))/data
mkdir -p test_${pop1}_${pop2}
cp ind_${pop1}.txt ind_${pop2}.txt test_${pop1}_${pop2}/
cd test_${pop1}_${pop2}/
cat ind_${pop1}.txt ind_${pop2}.txt > ind_${pop1}_${pop2}.txt
#计算SNP的fst
vcftools --gzvcf ../${fn}.vcf.gz --weir-fst-pop ind_${pop1}.txt --weir-fst-pop ind_${pop2}.txt --out ${pop1}_${pop2}
#筛选aim
cp ../../code/aim_select.py .
python aim_select.py --fst ${pop1}_${pop2}.weir.fst --ld ../${fn}_r2_morethan${r2_threshold}.ld --out aim_${pop1}_${pop2}.txt --fst-threshold ${fst_threshold} --ld-threshold ${r2_threshold}
#获取aim的SNP ID
awk 'NR==FNR {pos[$1":"$2]; next} ($1":"$4) in pos {print $2}' aim_${pop1}_${pop2}.txt ../${fn}.map > aim_${pop1}_${pop2}.snpid
paste aim_${pop1}_${pop2}.snpid aim_${pop1}_${pop2}.txt > tmp
rm aim_${pop1}_${pop2}.txt; mv tmp aim_${pop1}_${pop2}.txt

#提取个体 去死吧
awk '
    # 读取ind_${pop1}_${pop2}.txt文件，创建查找字典
NR==FNR {
    # 移除可能的换行符并存储
    gsub(/\r$/, "", $0)
    targets[$0] = 1
    next
}
    # 处理${fn}.fam文件
{
    # 构建"第一列_第二列"格式
    key = $2
    if (key in targets) {
        print $1 "\t" $2
    }
}' ind_${pop1}_${pop2}.txt ../${fn}.fam > ind_${pop1}_${pop2}.keep
plink --file ../${fn} --keep ind_${pop1}_${pop2}.keep --extract aim_${pop1}_${pop2}.snpid --recode --make-bed --out aim_${pop1}_${pop2}
#提取map的第2、4列 染色体名称、物理位置
cat aim_${pop1}_${pop2}.map | awk '{print $2"\t"$4}' > aim_${pop1}_${pop2}.info
cp ../../code/{haploview.sh,convert_plink2eigen.sh} .
#筛选tagSNP
> haploview.log
bash haploview.sh aim_${pop1}_${pop2} ${haploview_dir} 2>&1 | tee -a haploview.log

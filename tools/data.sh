#!/bin/bash

#cat /mnt/hdd1/HaodongChen/merged_new_guizhou51/smartpca_origin/merged_new_guizhou51_qc.fam | awk '{print $1}' | sort | uniq -c | sort -rnk1,1 | awk '{print $2}' > pop
#手动去除Yi Tujia Zhuang Mongol

#为了提取人群 根据PCA手动剔除离群个体
#grep -f pop /mnt/hdd1/HaodongChen/merged_new_guizhou51/smartpca_origin/merged_new_guizhou51_qc.fam | awk '{print $1"\t"$2}' > pop_sample_for_aim.txt

#plink --bfile /mnt/hdd1/HaodongChen/merged_new_guizhou51/smartpca_origin/merged_new_guizhou51_qc --keep pop_sample_for_aim.txt --make-bed --out merged_new

#按地区合并人群标签
#bash rename_tag.sh

#sh convert_eigen2plink.sh modern_EA_HO


# 质控
# qc
#plink --bfile modern_EA_HO --geno 0.1 --mind 0.1 --hwe 1e-3 --make-bed --out modern_EA_HO_qc
# 亲缘

# kingship
#king -b modern_EA_HO_qc.bed --kinship --degree 3
# 三级内亲缘关系对
#cat king.kin | awk '$9>0.0442' > relative_within3
##need to be removed individuals
##king --unrelate不行

#'''
#Atayal.HO	NA13606.HO
#Bonan.HO	BAO03.HO
#.....
#'''

# 共去除9个 实现：coze得到remove_ind.txt
#grep -f remove_ind.txt modern_EA_HO_qc.fam | awk '{print $1"\t"$2}' > remove_pop_ind.txt
#plink --bfile modern_EA_HO_qc --remove remove_pop_ind.txt --make-bed --out modern_EA_qc_kin3

# 改snp名字
#cat modern_EA_qc_kin3.bim | awk '$2=$1"_"$4{print}' OFS='\t' > tmp
#rm modern_EA_qc_kin3.bim ;mv tmp modern_EA_qc_kin3.bim

#合并标签后的 pop sample对应
cat modern_EA_qc_kin3.fam | awk '{print $1"\t"$2}' > pop_sample_for_aim_renametag.txt
bash pop_sample_file.sh

plink --bfile modern_EA_qc_kin3 --recode vcf-iid --out modern_EA_qc_kin3

bgzip modern_EA_qc_kin3.vcf
bcftools index modern_EA_qc_kin3.vcf.gz

bash convert_plink2eigen.sh modern_EA_qc_kin3

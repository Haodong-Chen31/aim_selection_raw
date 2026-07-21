#!/bin/sh

vcf_file=$1	#上游处理好的vcf文件
pop=$2		#对应的poplist文件
pop_sample=$3	#pop sample相对应的文件

fn=$(basename ${vcf_file} .alldone.vcf.gz)
code_dir=$(pwd)
work_dir=$(dirname ${code_dir})

cd ${work_dir}
##############################step1~5只是为了大致看人群结构，并确定用于位点筛选的人群和个体
##############################实际分成三个部分：
##############################第一step1~3——质控数据（包含geno,mind,LD），保证数据平衡性，初步PCA查看outlier
##############################第二step4~5——去除掉明显outlier和样本量小于5的群体，以减少其对后面AIM筛选的影响，剩余样本再次PCA观察，这些个体被归为core sample
##############################第三step6从只geno,mind质控后的数据中提取core sample，作为后续AIM筛选的样本
# step1: QC
sh ${code_dir}/step1_qc.sh ${vcf_file} ${fn} ${pop_sample} ${code_dir}

# step2: 对于样本量过大的人群，随机选择20个
#python ${code_dir}/sample_balance.py ${fn}_qc_pruned.fam outliers.txt balanced_samples.txt
#plink --bfile ${fn}_qc_pruned --keep balanced_samples.txt --make-bed --out ${fn}_balanced_subset
#sh ${code_dir}/convert_plink2eigen.sh ${fn}_balanced_subset

# step3: PCA看大致人群结构 并找出过滤掉outlier
sh ${code_dir}/smartpca.sh ${fn}_qc_pruned ${pop} ${code_dir}

# 查看outlier
python ${code_dir}/pca_plot_scatter_html_v2.py -p plot_${fn}_qc_pruned.txt -l ${pop}
mv scatter_plot.html smartpca_${fn}_qc_pruned.html

exit
# step4: 过滤outlier以及个体数小于5的人群
# 查看pca结果，编辑需要过滤掉的个体FID IID 存放到outlier_2nd_turn.txt
cd ${work_dir}
#提取上一步用于计算PC的样本
cat plot_${fn}_balanced_subset.txt | awk '{print $1,$2}' > keep_2nd_turn.txt
plink --bfile ${fn}_balanced_subset --keep keep_2nd_turn.txt --make-bed --out tmp_${fn}_balanced_subset
# 去除outlier_2nd_turn.txt，保留上一步用于计算的个体
plink --bfile tmp_${fn}_balanced_subset --remove outlier_2nd_turn.txt --make-bed --out ${fn}_balanced_subset_rm_outlier
rm tmp_${fn}_balanced_subset.*

# 统计样本量大于等于5的人群 FID IID sample_mt5.txt
python ${code_dir}/keep_sample_mt5.py ${fn}_balanced_subset_rm_outlier.fam ${fn}_balanced_subset_rm_outlier_sample_mt5.txt
plink --bfile ${fn}_balanced_subset_rm_outlier --keep ${fn}_balanced_subset_rm_outlier_sample_mt5.txt --make-bed --out ${fn}_balanced_subset_rm_outlier_sample_mt5
sh ${code_dir}/convert_plink2eigen.sh ${fn}_balanced_subset_rm_outlier_sample_mt5

# step5: 再跑一次PCA
# 需要注意——重新编辑人群列表
sh ${code_dir}/smartpca.sh ${fn}_balanced_subset_rm_outlier_sample_mt5 ${pop}
python ${code_dir}/pca_plot_scatter_html_v2.py -p plot_${fn}_balanced_subset_rm_outlier_sample_mt5.txt -l ${pop}
mv scatter_plot.html smartpca_${fn}_balanced_subset_rm_outlier_sample_mt5.html

#结果好，就确定core sample
cat plot_${fn}_balanced_subset_rm_outlier_sample_mt5.txt | awk '{print $1,$2}' > core_samples.txt
#exit
# step6: 
plink --bfile ${fn}_qc --keep core_samples.txt --make-bed --out core_${fn}


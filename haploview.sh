#!/bin/bash
#Haploview 4.2挑选tagSNP

prefix_file=$1
haploview_dir=$2

#JAVA8=/opt/java/zulu8.23.0.3-jdk8.0.144-linux_x64/bin/java
#export prefix_file
#export haploview_dir

#export JAVA_HOME=/opt/java/zulu8.23.0.3-jdk8.0.144-linux_x64
#export PATH=$JAVA_HOME/bin:$PATH
# 数据包含多个染色体，分开处理
#for chr in {1..22}; do
#    # 提取单个染色体的数据
#    plink --file ${prefix_file} --chr ${chr} --threads 16 --recode --out chr${chr} > /dev/null 2>&1
#    cat chr${chr}.map | awk '{print $2"\t"$4}' > chr${chr}.info
#    echo "Haploview-chr${chr}运行开始，输入SNP数量：$(wc -l < chr${chr}.info)"
#    # 运行Haploview
#    java -Xmx4g -jar ${haploview_dir}/Haploview.jar -n -pedfile chr${chr}.ped -info chr${chr}.info -pairwiseTagging -skipcheck -dprime -tagrsqcutoff 0.8 -minMAF 0.05 -out chr${chr}_tagSNPs > /dev/null 2>&1
#    echo "Haploview-chr${chr}运行完成，tagSNP数量：$(wc -l < chr${chr}_tagSNPs.TESTS)"
#done

#parallel -j 5 '
#    chr={}
#
#    plink --file ${prefix_file} --chr ${chr} --threads 4 \
#          --recode --out chr${chr} > /dev/null 2>&1
#
#    awk "{print \$2\"\t\"\$4}" chr${chr}.map > chr${chr}.info
#
#    echo "Haploview-chr${chr} 运行开始，输入SNP数量：$(wc -l < chr${chr}.info)"
#
#    java -Xmx4g -jar ${haploview_dir}/Haploview.jar \
#         -n -pedfile chr${chr}.ped \
#         -info chr${chr}.info \
#         -pairwiseTagging -skipcheck -dprime \
#         -tagrsqcutoff 0.8 -minMAF 0.05 \
#         -out chr${chr}_tagSNPs
#
#    echo "Haploview-chr${chr} 运行完成，tagSNP数量：$(wc -l < chr${chr}_tagSNPs.TESTS)"
#' ::: {1..22}

# LD过滤
plink --bfile ${prefix_file} --indep-pairwise 500 25 0.1 --out ${prefix_file}
plink --bfile ${prefix_file} --extract ${prefix_file}.prune.in --make-bed --recode --out tagSNP_${prefix_file#aim_}

cat tagSNP_${prefix_file#aim_}.bim | awk '{print $2}'> tagSNP_${prefix_file#aim_}.snpid

echo "LD过滤完成，tagSNP数量：$(wc -l < tagSNP_${prefix_file#aim_}.snpid)"
#rm -rf chr*
awk 'NR==FNR{a[$1]; next} ($1 in a)' tagSNP_${prefix_file#aim_}.snpid ${prefix_file}.txt > tagSNP_${prefix_file#aim_}.txt
#plink --file ${prefix_file} --extract tagSNP_${prefix_file#aim_}.snpid --recode --make-bed --out tagSNP_${prefix_file#aim_} > /dev/null 2>&1

#top

#质控
echo "质控（--geno 0.1 --mind 0.1 --hwe 1e-3）..."
plink --bfile ${prefix_file} --geno 0.1 --mind 0.1 --hwe 1e-3 --make-bed --out ${prefix_file}_qc > /dev/null 2>&1
plink --bfile tagSNP_${prefix_file#aim_} --geno 0.1 --mind 0.1 --hwe 1e-3 --make-bed --out tagSNP_${prefix_file#aim_}_qc > /dev/null 2>&1
echo "质控完成"

#bed/bim/fam -> geno/snp/ind
echo "bed/bim/fam -> geno/snp/ind..."
bash convert_plink2eigen.sh ${prefix_file} > /dev/null 2>&1
bash convert_plink2eigen.sh ${prefix_file}_qc > /dev/null 2>&1
bash convert_plink2eigen.sh tagSNP_${prefix_file#aim_} > /dev/null 2>&1
bash convert_plink2eigen.sh tagSNP_${prefix_file#aim_}_qc > /dev/null 2>&1
echo "格式转换完成"

echo "==tagSNP筛选完毕！=="
echo "AIM数量：		$(wc -l < ${prefix_file}.map)，文件：${prefix_file}.{snpid,bed/bim/fam,ped/map,geno/snp/ind}"
echo "质控后AIM数量：		$(wc -l < ${prefix_file}_qc.bim)，文件：${prefix_file}_qc.{bed/bim/fam,geno/snp/ind}"
echo "tagSNP数量：		$(wc -l < tagSNP_${prefix_file#aim_}.map)，文件：tagSNP_${prefix_file#aim_}.{snpid,bed/bim/fam,ped/map,geno/snp/ind}"
echo "质控后tagSNP数量：	$(wc -l < tagSNP_${prefix_file#aim_}_qc.bim)，文件：tagSNP_${prefix_file#aim_}_qc.{bed/bim/fam,geno/snp/ind}"

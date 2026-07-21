#!/bin/sh

fn=merged_1kg_sgdp.qc

mkdir -p $(dirname $(pwd))/log
cp sample_ind.sh ../data
cd $(dirname $(pwd))/data/
#获得每个人群包含的样本名
sh sample_ind.sh
cat ind_* > all_EA_samples.txt
ls | egrep ind_ | sed 's/^ind_//g ; s/.txt$//g' > pop
#提东亚人群
#bcftools view -S all_EA_samples.txt merged_autosome_qc.vcf.gz -Oz -o ${fn}.vcf.gz
#sort
bcftools sort ${fn}.vcf.gz -Oz -o ${fn}.sorted.vcf.gz
bcftools index ${fn}.sorted.vcf.gz

#改SNP的ID
zcat ${fn}.vcf.gz | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$3 = $1 "_" $2; print}' | bgzip > ${fn}_snp_renamed.vcf.gz
rm ${fn}.vcf.gz ; mv ${fn}_snp_renamed.vcf.gz ${fn}.vcf.gz
bcftools index ${fn}.vcf.gz
plink --vcf ${fn}.vcf.gz --make-bed --out ${fn}
# sort
#bcftools sort ${fn}.vcf.gz -Oz -o ${fn}.sorted.vcf.gz
# 修改fam
> ind2pop.txt
for f in $(ls ind_*.txt);do pop=$(echo ${f} | sed 's/^ind_//g;s/.txt$//g') ; awk -v pop=${pop} '{print $1"\t"pop}' ${f} >> ind2pop.txt;done
awk 'BEGIN{OFS="\t"}
     NR==FNR {pop[$1]=$2; next}
     {
       if ($2 in pop) $1=pop[$2];
       else $1="UNKNOWN";
       print
     }' ind2pop.txt ${fn}.fam > tmp.fam
rm ${fn}.fam;mv tmp.fam ${fn}.fam
plink --bfile ${fn} --recode --out ${fn}
cat ${fn}.ped | awk '$6=1{print}' > tmp.ped
rm ${fn}.ped ; mv tmp.ped ${fn}.ped

cd /mnt/hdd7/HaodongChen/aim_260106_liftOver/code/
sh step1_aim_tagSNP_select.parl.sh ${fn} 500 0.05 0.1 ../data/pop 32
sh step2_test4.parl.sh ${fn} | tee step2.log
sh combined_fst.sh ../data/pop ${fn} > ../log/combined_fst.log
sh snp_statistics.sh

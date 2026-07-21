#!/bin/sh

fn_input=$1

cd $(dirname $(pwd))/data
cp ../code/convert_plink2eigen.sh .
for m in aim tagSNP;do
    for n in intersection union;do
	echo "===${m}_${n}数据集==="
        plink --bfile ${fn_input} --extract ${m}_${n}.snpid --make-bed --out ${m}_${n} > /dev/null 2>&1
	echo "质控（--geno 0.1 --mind 0.1 --hwe 1e-3）..."
        plink --bfile ${m}_${n} --geno 0.1 --mind 0.1 --hwe 1e-3 --make-bed --out ${m}_${n}_qc > /dev/null 2>&1
	echo "格式转换 {bed,bim,fam} -> {geno,snp,ind}..."
        sh convert_plink2eigen.sh ${m}_${n} > /dev/null 2>&1
        sh convert_plink2eigen.sh ${m}_${n}_qc > /dev/null 2>&1
	echo "${m}_${n}数据集完成！"
    done
done



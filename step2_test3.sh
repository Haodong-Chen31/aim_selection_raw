#!/bin/sh

mapfile -t pops < $1
mean_threshold=$2
cv_threshold=$3

fst_files=()
pop_pairs=()

cd $(dirname $(pwd))/data
# 收集文件和对应的种群对名称
for ((i=0; i<${#pops[@]}; i++)); do
    for ((j=i+1; j<${#pops[@]}; j++)); do
        pop_pair="${pops[i]}_${pops[j]}"
        fst_file="./test_${pop_pair}/${pop_pair}.weir.fst"

        if [[ -f ${fst_file} ]]; then
            fst_files+=("${fst_file}")
            pop_pairs+=("${pop_pair}")
            echo "找到文件: ${fst_file}"
        else
            echo "错误: 文件不存在 - ${fst_file}"
        fi
    done
done > fstfile_exists_or_not.txt

# 检查bim文件是否存在
bim_file="modern_EA_qc_kin3.bim"
if [[ ! -f ${bim_file} ]]; then
    echo "错误: bim文件不存在 - ${bim_file}"
    exit 1
fi

# 合并文件
if [[ ${#fst_files[@]} -gt 0 ]]; then
    {
        # 输出标题行
        echo -n "CHROM  POS     SNPID   ${pop_pairs[0]}"
        for ((i=1; i<${#pop_pairs[@]}; i++)); do
            echo -n "   ${pop_pairs[i]}"
        done
        echo ""

        # 输出数据行，添加SNP ID
        paste "${fst_files[@]}" | awk -v count=${#fst_files[@]} -v bim_file="${bim_file}" '
        BEGIN {
            # 读取bim文件，建立位置到SNP ID的映射
            while ((getline < bim_file) > 0) {
                # 使用染色体和位置作为键，SNP ID作为值
                key = $1 "\t" $4
                snp_id[key] = $2
            }
            close(bim_file)
        }
        NR>1 {
            chrom = $1
            pos = $2
            # 查找对应的SNP ID
            key = chrom "\t" pos
            if (key in snp_id) {
                snpid = snp_id[key]
            } else {
                snpid = "NA"
            }
            printf "%s\t%s\t%s\t%s", chrom, pos, snpid, $3
            for(i=2; i<=count; i++) {
                printf "\t%s", $(i*3)
            }
            printf "\n"
        }'
    } > combined_fst.txt

    echo "合并完成！共处理 ${#fst_files[@]} 个种群对的FST数据"
    echo "输出文件: combined_fst.txt"
    echo "列信息: CHROM, POS, SNPID, ${pop_pairs[0]}, ${pop_pairs[1]}, ..."
else
    echo "错误: 没有找到任何FST文件"
    exit 1
fi

#上面代码不需要每次都跑
#根据mean_threshold和cv_threshold选SNP
Rscript ../code/SNPselect_method1.r combined_fst.txt ${mean_threshold} ${cv_threshold}
cp ../code/convert_plink2eigen.sh .
plink --bfile modern_EA_qc_kin3 --extract combined_fst_candidate_snps_mean${mean_threshold}_cv${cv_threshold}.snpid --make-bed --out combined_fst
plink --bfile combined_fst --geno 0.1 --mind 0.1 --hwe 1e-3 --make-bed --out combined_fst_qc
sh convert_plink2eigen.sh combined_fst_qc

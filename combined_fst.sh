#!/usr/bin/env bash

mapfile -t pops < $1
prefix_fn=$2
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
bim_file="${prefix_fn}.bim"
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

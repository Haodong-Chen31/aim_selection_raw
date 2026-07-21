#!/bin/sh

mapfile -t pops < $1
aim_pop1_pop2=""
tagSNP_pop1_pop2=""
missing_files=""
all_files_exist=true

cd $(dirname $(pwd))/data
for ((i=0; i<${#pops[@]}; i++)); do
    for ((j=i+1; j<${#pops[@]}; j++)); do
        # 定义当前种群对的文件名
        aim_file="./test_${pops[i]}_${pops[j]}/aim_${pops[i]}_${pops[j]}.snpid"
        tag_file="./test_${pops[i]}_${pops[j]}/tagSNP_${pops[i]}_${pops[j]}.snpid"
        
        # 检查是否每组都生成了对应的文件
        if [[ -f ${aim_file} && -f ${tag_file} ]]; then
            # 将文件名添加到数组中
            aim_pop1_pop2="${aim_pop1_pop2} ${aim_file}"
            tagSNP_pop1_pop2="${tagSNP_pop1_pop2} ${tag_file}"
            echo "找到文件: ${aim_file} 和 ${tag_file}"
        else
            echo "错误: 文件不存在 - ${aim_file} 或 ${tag_file}"
            all_files_exist=false
            missing_files="${missing_files} ${aim_file} ${tag_file}"
        fi
    done
done > file_exists_or_not.txt

# 如果有任何文件缺失，则退出不执行后续代码
if [[ ${all_files_exist} != "true" ]]; then
    echo "错误: 以下文件缺失: ${missing_files}"
    echo "终止执行，不进行SNP交集和并集计算"
    exit 1
fi

# 所有文件都存在，执行SNP的交集和并集计算
echo "所有文件都存在，开始计算SNP的交集和并集..."
echo "找到 $(echo ${aim_pop1_pop2} | wc -w) 个AIM文件和 $(echo ${tagSNP_pop1_pop2} | wc -w) 个tagSNP文件"
#交集
sort ${aim_pop1_pop2} | uniq -d > aim_intersection.snpid
sort ${tagSNP_pop1_pop2} | uniq -d > tagSNP_intersection.snpid
#并集
sort ${aim_pop1_pop2} | uniq > aim_union.snpid
sort ${tagSNP_pop1_pop2} | uniq > tagSNP_union.snpid

echo "AIM数量（交集）：    $(wc -l < aim_intersection.snpid)"
echo "AIM数量（并集）：    $(wc -l < aim_union.snpid)"
echo "tagSNP数量（交集）： $(wc -l < tagSNP_intersection.snpid)"
echo "tagSNP数量（并集）： $(wc -l < tagSNP_union.snpid)"

#取SNP的交集和并集
#交集
#sort --parallel=32 -T ~/tmp -S 4G --mmap ${aim_pop1_pop2} | uniq -d > aim_intersection.snpid
#sort --parallel=32 -T ~/tmp -S 4G --mmap ${tagSNP_pop1_pop2} | uniq -d > tagSNP_intersection.snpid
#并集
#sort --parallel=32 -T ~/tmp -S 4G --mmap ${aim_pop1_pop2} | uniq > aim_union.snpid
#sort --parallel=32 -T ~/tmp -S 4G --mmap ${tagSNP_pop1_pop2} | uniq > tagSNP_union.snpid

#!/bin/bash

mapfile -t pops < $1
top_num=$2

dir=$(dirname $(pwd))
for ((i=0; i<${#pops[@]}; i++)); do
    for ((j=i+1; j<${#pops[@]}; j++)); do
        cd ${dir}/data/test_${pops[i]}_${pops[j]}
        paste aim_${pops[i]}_${pops[j]}.snpid aim_${pops[i]}_${pops[j]}.txt > tmp
        rm aim_${pops[i]}_${pops[j]}.txt; mv tmp aim_${pops[i]}_${pops[j]}.txt
        awk 'NR==FNR{a[$1]; next} ($1 in a)' tagSNP_${pops[i]}_${pops[j]}.snpid aim_${pops[i]}_${pops[j]}.txt > tagSNP_${pops[i]}_${pops[j]}.txt
        cat tagSNP_${pops[i]}_${pops[j]}.txt | sort -rk4,4 | head -${top_num} | awk '{print $1}' > tagSNP_${pops[i]}_${pops[j]}.top${top_num}.snpid
    done
done

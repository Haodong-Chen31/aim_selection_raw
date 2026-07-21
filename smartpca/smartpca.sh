#!/bin/sh

fn=$1
pop=$2
geno_dir=$3
code_dir=$4

#alias rmsp='sed "s/^\s*//g" | sed "s/[[:blank:]]\+/\t/g"'

rmsp() {
    sed 's/^[[:space:]]*//g' | sed 's/[[:blank:]]\+/\t/g'
}

# smartpca.par
for i in "smartpca_${fn}.par";do
    echo "genotypename: ${geno_dir}/${fn}.geno"
    echo "snpname:      ${geno_dir}/${fn}.snp"
    echo "indivname:    ${geno_dir}/${fn}.ind"
    echo "evecoutname:  smartpca_${fn}.evec"
    echo "evaloutname:  smartpca_${fn}.eval"
    echo "poplistname:  ${pop}"
    echo "lsqproject: YES"
    echo "numoutevec: 5"
    echo "altnormstyle: NO"
    echo "numoutlieriter: 0"
    echo "numthreads: 4"
done > smartpca_${fn}.par

# smartpca
smartpca -p smartpca_${fn}.par > smartpca_${fn}.log 2>&1

# Calculate PCs
lines=$(wc -l smartpca_${fn}.eval | rmsp | cut -f 1)
lines=$(expr ${lines} - 1 )
pc1=$(head -n 1 smartpca_${fn}.eval)
pc2=$(tail -n+2 smartpca_${fn}.eval | head -n 1)
echo -n "PC1 (" >  PCs_${fn}.txt ; echo "${pc1}/${lines}*100" | xargs -n 1 python ${code_dir}/bc.py >> PCs_${fn}.txt ; echo "%)" >> PCs_${fn}.txt
echo -n "PC2 (" >> PCs_${fn}.txt ; echo "${pc2}/${lines}*100" | xargs -n 1 python ${code_dir}/bc.py >> PCs_${fn}.txt ; echo "%)" >> PCs_${fn}.txt
cat smartpca_${fn}.log | grep "total number of snps killed in pass" | tail -n 1                                >> PCs_${fn}.txt

# Post-Processing
tail -n+2 smartpca_${fn}.evec | awk '{print $7,$1,$2,$3}' > plot_${fn}.txt
python ${code_dir}/pcaRploter_v3.py -p ${pop} -o smartpca_${fn}.r -i plot_${fn}.txt --pdf smartpca_${fn}.pdf --legend smartpca_${fn}_legend.pdf
Rscript smartpca_${fn}.r
#zip ${fn}.zip *.{py,pdf,r,sh,txt} smartpca.* *poplist*

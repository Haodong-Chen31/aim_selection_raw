#!/usr/bin/env bash

#usage: sh convert.sh [prefix of file]
DIR=$(pwd)
fn=$1

PARFILE=$DIR/parconvert_eigen2plink_${fn}
cat >$PARFILE<<EOF
genotypename:	${DIR}/${fn}.geno
snpname:	${DIR}/${fn}.snp
indivname:	${DIR}/${fn}.ind
outputformat:	PACKEDPED
genooutfilename:	${DIR}/${fn}.bed
snpoutfilename:	${DIR}/${fn}.bim
indoutfilename:	${DIR}/${fn}.fam
EOF

convertf -p $PARFILE

#.fam文件改对
cat ${DIR}/${fn}.ind > ${DIR}/1.txt
printf "\n" >> ${DIR}/1.txt;cat ${DIR}/${fn}.ind >> ${DIR}/1.txt
awk -F " " 'NR==FNR{pop[$1]=$3}NR>FNR{print pop[$2],substr($0,index($0,$2))}' ${DIR}/1.txt ${DIR}/${fn}.fam > ${DIR}/${fn}_1.fam
rm ${DIR}/${fn}.fam ${DIR}/1.txt
mv ${DIR}/${fn}_1.fam ${DIR}/${fn}.fam

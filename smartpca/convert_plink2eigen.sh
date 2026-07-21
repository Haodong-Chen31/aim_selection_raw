#!/bin/bash

#usage: sh convert.sh [prefix of file]
DIR=$(pwd)
fn=$1

PARFILE=$DIR/parconvert_plink2eigen_${fn}
cat<<EOF>$PARFILE
genotypename:	${DIR}/${fn}.bed
snpname:	${DIR}/${fn}.bim
indivname:	${DIR}/${fn}.fam
outputformat:	EIGENSTRAT
genooutfilename:	${DIR}/${fn}.geno
snpoutfilename:	${DIR}/${fn}.snp
indoutfilename:	${DIR}/${fn}.ind
EOF
convertf -p $PARFILE

####################################
#.ind文件改对
sed 's/[ ][ ]*/ /g' ${DIR}/${fn}.ind | sed 's/:/ /' | awk '{print $2" "$3" "$1}' > ${DIR}/${fn}_tmp.ind

rm ${DIR}/${fn}.ind;mv ${DIR}/${fn}_tmp.ind ${DIR}/${fn}.ind

####################################
#改.snp文件
cat ${fn}.snp | awk '{print $2"_"$4"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ${fn}_tmp.snp

rm ${fn}.snp; mv ${fn}_tmp.snp ${fn}.snp

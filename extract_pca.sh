#!/usr/bin/bash
#usage: sh extract_pca.sh [merged file name(path)] [poplist for pca] [modern pop]

DIR=$(pwd)
fn=$1
POP=$2
POP_M=${POP}

##extract
f=$(basename ${fn})
cat >parextract_${f}<<EOF
genotypename:	${fn}.geno
snpname:	${fn}.snp
indivname:	${fn}.ind
outputformat:	EIGENSTRAT
genooutfilename:	${DIR}/$(basename ${fn}).geno
snpoutfilename:	${DIR}/$(basename ${fn}).snp
indoutfilename:	${DIR}/$(basename ${fn}).ind
poplistname:	${DIR}/${POP}
hashcheck:	NO
strandcheck:	NO
allowdups:	YES
EOF
convertf -p parextract_${f}

##smartpca
cat >parsmartpca_${f}<<EOF
genotypename:	${DIR}/$(basename ${fn}).geno
snpname:	${DIR}/$(basename ${fn}).snp
indivname:	${DIR}/$(basename ${fn}).ind
evecoutname:	${DIR}/$(basename ${fn}).evec
evaloutname:	${DIR}/$(basename ${fn}).eval
poplistname:	${DIR}/${POP_M}
lsqproject:	YES
numoutevec:	10
numoutlieriter:	0
altnormstyle:	NO
EOF
smartpca -p parsmartpca_${f}

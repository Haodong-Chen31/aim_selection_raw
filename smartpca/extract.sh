#!/bin/bash

#usage: sh ~/code/convert_merge/extract.sh [prefix of filename (path)] [target directory] [prefix of target filename] [poplist]

DIR=$(pwd)
fn=$1
fn_out=$2
pop=$3

cat<<EOF>parextract
genotypename:	${fn}.geno
snpname:	${fn}.snp
indivname:	${fn}.ind
outputformat:	EIGENSTRAT
genooutfilename:	${DIR}/${fn_out}.geno
snpoutfilename:	${DIR}/${fn_out}.snp
indoutfilename:	${DIR}/${fn_out}.ind
poplistname:	${DIR}/${pop}
hashcheck:	NO
strandcheck:	NO
allowdups:	YES
EOF

convertf -p parextract

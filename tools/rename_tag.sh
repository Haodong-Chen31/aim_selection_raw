#!/bin/bash

for f in $(ls *.fam);do
	sed -i -e 's/\(Han_Fujian\|Han_Hubei\|Han_Jiangsu\|Han_Sichuan\|Han_Zhejiang\)/Han_South/g' \
       -e 's/\(Han_Henan\|Han_Shandong\|Han_Shanxi\)/Han_North/g' \
       -e 's/\(Tibetan_Chamdo\|Tibetan_Lhasa\|Tibetan_Nagqu\|Tibetan_Shannan\)/Tibetan_Xizang/g' \
       -e 's/\(Dai.HO\|Mulam.HO\|Maonan.HO\|Li.HO\|Nung\)/Tai_Kadai/g' \
       -e 's/\(She.HO\|Miao.HO\)/Hmong_Mien/g' \
       -e 's/\(Ami.HO\|Atayal.HO\)/Austronesian/g' \
       -e 's/Mongol.HO/Mongolic/g' \
       -e 's/Ulchi.HO/Tungusic/g' \
       ${f}
done

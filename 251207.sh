#!/bin/sh

sh step2_test4.parl.sh
wait
sh snp_statistics.sh
wait
cd ../smartpca_more_extract
sh pca.parl.sh

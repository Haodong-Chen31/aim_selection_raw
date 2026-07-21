#!/bin/bash

for n in $(seq 0 100 900);do
	Rscript nb_mcc_functions.R $((n+1)) $((n+100))
done

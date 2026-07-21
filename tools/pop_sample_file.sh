#!/bin/sh

awk '{print $2 > "ind_" $1 ".txt"}' pop_sample_for_aim.txt

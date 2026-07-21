# aim_selection_raw
A pipeline for ancestry informative SNP panel based on FST and LD, just for personal use.

# 用到的软件, R包, python库等
```
#PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)
#VCFtools (0.1.16)
#python 3.13.5
#import pandas as pd
#import argparse
#from collections import defaultdict
#import sys
#openjdk 21.0.6 2025-01-21
#OpenJDK Runtime Environment JBR-21.0.6+9-895.97-nomod (build 21.0.6+9-b895.97)
#OpenJDK 64-Bit Server VM JBR-21.0.6+9-895.97-nomod (build 21.0.6+9-b895.97, mixed mode, sharing)
#Haploview.jar
#R version 4.3.1 (2023-06-16)
#library(qqman)
library(caret)
library(naivebayes)
library(ggplot2)
library(reshape2)
library(genio)
library(patchwork)
library(ggthemes)
library(foreach)
library(doParallel)
```

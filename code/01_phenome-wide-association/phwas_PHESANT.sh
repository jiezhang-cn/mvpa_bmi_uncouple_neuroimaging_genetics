#!/bin/bash

nohup Rscript /home/user2/jzhang_data/Dissertation/code/PHESANT/WAS/phenomeScan.r \
  --phenofile="/home/user2/jzhang_data/Dissertation/data/ukb_data_for_PHESANT.csv" \
  --traitofinterestfile="/home/user2/jzhang_data/Dissertation/data/PHESCANT_pa_obesity_pheno.csv" \
  --variablelistfile="/home/user2/jzhang_data/Dissertation/code/PHESANT/variable-info/outcome-info.tsv" \
  --datacodingfile="/home/user2/jzhang_data/Dissertation/code/PHESANT/variable-info/data-coding-ordinal-info.txt" \
  --traitofinterest="high_mvpa_bmi_uncouple" \
  --resDir="/home/user2/jzhang_data/Dissertation/results/high_pa_obesity_uncouple_results" \
  --confounderfile="/home/user2/jzhang_data/Dissertation/data/PHESCANT_pa_obesity_covariate.csv" \
  --numParts=20
  > /home/user2/jzhang_data/Dissertation/code/PHESANT_high_mvpa.log 2>&1 &
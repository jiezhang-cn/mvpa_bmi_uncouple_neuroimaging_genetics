#!/bin/bash

GCTB="/home/user2/jzhang_data/code/gctb/gctb"
ANNOT="/home/user2/jzhang_data/code/gctb_ref/annot_bolt_clean.txt"
GENEMAP="/home/user2/jzhang_data/code/gctb_ref/gene_map_hg38_hg19.txt"

# ---- high_mvpa ----
LDM_EIGEN1="/test/jzhang_data_depository/gctb_finemapping_results/matched_ldm"
GWAS1="/home/Disk/Data2/imputation.bgen/PA_obesity_uncouple_GWAS/high_mvpa_bmi_uncouple.ma"
OUT1="/test/jzhang_data_depository/gctb_finemapping_results/high_mvpa_gwfm"

${GCTB} --gwfm RC \
  --ldm-eigen ${LDM_EIGEN1} \
  --gwas-summary ${GWAS1} \
  --annot ${ANNOT} \
  --gene-map ${GENEMAP} \
  --thread 56 \
  --out ${OUT1} \
  > ${OUT1}.log 2>&1

echo "high_mvpa done, exit code: $?"

# ---- low_mvpa ----
LDM_EIGEN2="/test/jzhang_data_depository/gctb_finemapping_results2/matched_ldm"
GWAS2="/home/Disk/Data2/imputation.bgen/PA_obesity_uncouple_GWAS/low_mvpa_bmi_uncouple.ma"
OUT2="/test/jzhang_data_depository/gctb_finemapping_results2/low_mvpa_gwfm"

${GCTB} --gwfm RC \
  --ldm-eigen ${LDM_EIGEN2} \
  --gwas-summary ${GWAS2} \
  --annot ${ANNOT} \
  --gene-map ${GENEMAP} \
  --thread 56 \
  --out ${OUT2} \
  > ${OUT2}.log 2>&1

echo "low_mvpa done, exit code: $?"
echo "All done"

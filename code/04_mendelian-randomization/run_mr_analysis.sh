#!/bin/bash
# ============================================================================
#  Bi-directional MR + CAUSE pipeline
#  Decoupled MVPA phenotypes (BMI-uncoupled) <-> Psychiatric disorders
# ============================================================================
set -uo pipefail

BASE=/home/user2/jzhang_data
OUT=${BASE}/results/mr_results
mkdir -p ${OUT}/{exposure_dat,outcome_dat,harmonised,results,sensitivity,plots,logs,cache}

MVPA_HIGH=/home/Disk/Data2/imputation.bgen/PA_obesity_uncouple_GWAS/high_mvpa_bmi_uncouple_autosomes_bgen_stats.gz
MVPA_LOW=/home/Disk/Data2/imputation.bgen/PA_obesity_uncouple_GWAS/low_mvpa_bmi_uncouple_autosomes_bgen_stats.gz

GWAS_DIR=${BASE}/data/psychiatric_trait_gwas
CD_DIR=${GWAS_DIR}/Cross_disorder
PREP_DIR=${BASE}/results/ldsc_full/prepared      # cleaned files from the LDSC stage, pre-munge
REF_BFILE=${BASE}/code/reference_data/g1000_eur
PLINK_BIN=${BASE}/code/plink

echo "==================================================="
echo "  Bi-directional MR pipeline: $(date)"
echo "==================================================="

Rscript ${BASE}/code/mr_pipeline.R \
    --out          ${OUT} \
    --mvpa_high    ${MVPA_HIGH} \
    --mvpa_low     ${MVPA_LOW} \
    --gwas_dir     ${GWAS_DIR} \
    --cd_dir       ${CD_DIR} \
    --prep_dir     ${PREP_DIR} \
    --bfile        ${REF_BFILE} \
    --plink_bin    ${PLINK_BIN} \
    2>&1 | tee ${OUT}/logs/mr_pipeline.log

echo "==================================================="
echo "  Plotting summary"
echo "==================================================="

Rscript ${BASE}/code/plot_mr_summary.R --out ${OUT} \
    2>&1 | tee ${OUT}/logs/plot_summary.log

echo "==================================================="
echo "  Done: $(date)"
echo "==================================================="

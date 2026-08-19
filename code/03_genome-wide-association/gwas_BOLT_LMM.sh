#!/bin/bash

BOLT_CMD="/home/Disk/Data2/imputation.bgen/BOLT-LMM_v2.5/bolt"

LDSCORE_FILE="/home/Disk/Data2/imputation.bgen/BOLT-LMM_v2.5/tables/LDSCORE.1000G_EUR.tab.gz" 
GENETIC_MAP="/home/Disk/Data2/imputation.bgen/BOLT-LMM_v2.5/tables/genetic_map_hg19_withX.txt.gz"

IMPUTE_DIR="/home/Disk/Data2/imputation.bgen"
PHENO_FILE="${IMPUTE_DIR}/PA_obesity_uncouple_GWAS/Pheno_mvpa_bmi_uncouple"
COVAR_FILE="${IMPUTE_DIR}/PA_obesity_uncouple_GWAS/Covariate_mvpa_bmi_uncouple"
OUTPUT_DIR="${IMPUTE_DIR}/PA_obesity_uncouple_GWAS"

SAMPLE_FILE="${IMPUTE_DIR}/ukb22828_c1_b0_v3_s487283.sample"

PLINK_BED="/home/Disk/Data2/GWAS_analysis/ukb22418_c{1:22}_b0_v2.bed"
PLINK_BIM="/home/Disk/Data2/GWAS_analysis/ukb_snp_chr{1:22}_v2.bim"
PLINK_FAM="/home/Disk/Data2/GWAS_analysis/ukb22418_c1_b0_v2_s488251.fam" 

PHENOTYPES=("high_mvpa_bmi_uncouple" "low_mvpa_bmi_uncouple")

mkdir -p "$OUTPUT_DIR"

for PHENO in "${PHENOTYPES[@]}"; do
    echo "=================================================="
    echo "Processing phenotype: $PHENO"
    echo "=================================================="
    $BOLT_CMD \
        --bed="$PLINK_BED" \
        --bim="$PLINK_BIM" \
        --fam="$PLINK_FAM" \
        --remove=/home/Disk/Data2/imputation.bgen/PA_obesity_uncouple_GWAS/bolt.in_plink_but_not_imputed.FID_IID.968.txt \
        --LDscoresFile="$LDSCORE_FILE" \
        --geneticMapFile="$GENETIC_MAP" \
        --phenoFile="$PHENO_FILE" \
        --phenoCol="$PHENO" \
        --covarFile="$COVAR_FILE" \
        --covarCol=sex \
        --qCovarCol=age \
        --qCovarCol=age_2 \
        --qCovarCol=age_sex \
        --qCovarCol=pc{1:20} \
        --lmm \
        --lmmForceNonInf \
        --numThreads=24 \
        --statsFile="${OUTPUT_DIR}/${PHENO}_autosomes_stats.gz" \
        --bgenFile="${IMPUTE_DIR}/ukb_imp_chr{1:22}_v3.bgen" \
        --bgenMinMAF=1e-3  \
        --bgenMinINFO=0.3 \
        --sampleFile="$SAMPLE_FILE" \
        --statsFileBgenSnps="${OUTPUT_DIR}/${PHENO}_autosomes_bgen_stats.gz" \
        --verboseStats

done
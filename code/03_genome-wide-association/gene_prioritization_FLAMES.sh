#!/bin/bash
set -euo pipefail

# ========= 路径配置 =========
FLAMES_DIR="/home/user2/jzhang_data/Dissertation/code/FLAMES"
ANNOT_DIR="/home/user2/jzhang_data/Dissertation/code/FLAMES/annotation_file/Annotation_data"

PROJECT_DIR="/home/user2/jzhang_data/Dissertation/results/high_pa_obesity_uncouple_results"
POPS_PREDS="${PROJECT_DIR}/pops_results/high_mvpa_pops.preds"
MAGMA_OUT="${PROJECT_DIR}/magma_results/step2_gene_high_mvpa.genes.out"
MAGMA_TISSUE="${PROJECT_DIR}/magma_results/high_mvpa_tissue.gsa.out"
INDEX_FILE="${PROJECT_DIR}/flames_input/indexfile.txt"
FLAMES_OUT_DIR="${PROJECT_DIR}/flames_output"

# ========= VEP & CADD 本地配置 =========
REF_DIR="/home/user2/jzhang_data/Dissertation/code/flames_ref_data"
VEP_CACHE="${REF_DIR}/vep_cache"
CADD_FILE="${REF_DIR}/whole_genome_SNVs.tsv.gz"
TABIX_CMD="${HOME}/miniforge3/envs/FLAMES/bin/tabix"
VEP_CMD="${HOME}/miniforge3/envs/FLAMES/bin/vep"
# ============================

mkdir -p "$FLAMES_OUT_DIR"

# Step 1: Annotate（使用本地 VEP 和 CADD）
echo "[$(date)] Running FLAMES annotate..."
python "${FLAMES_DIR}/FLAMES.py" annotate \
    -a  "$ANNOT_DIR" \
    -p  "$POPS_PREDS" \
    -m  "$MAGMA_OUT" \
    -mt "$MAGMA_TISSUE" \
    -id "$INDEX_FILE" \
    -cv  "$VEP_CMD" \
    -vc  "$VEP_CACHE" \
    -t  "$TABIX_CMD" \
    -cf "$CADD_FILE"

# Step 2: Scoring
echo "[$(date)] Running FLAMES scoring..."
python "${FLAMES_DIR}/FLAMES.py" FLAMES \
    -id "$INDEX_FILE" \
    -o  "$FLAMES_OUT_DIR"

echo "[$(date)] Done. Output: $FLAMES_OUT_DIR"



#!/bin/bash
set -euo pipefail

# ========= 路径配置 =========
FLAMES_DIR="/home/user2/jzhang_data/Dissertation/code/FLAMES"
ANNOT_DIR="/home/user2/jzhang_data/Dissertation/code/FLAMES/annotation_file/Annotation_data"

PROJECT_DIR="/home/user2/jzhang_data/Dissertation/results/low_pa_obesity_uncouple_results"
POPS_PREDS="${PROJECT_DIR}/pops_results/low_mvpa_pops.preds"
MAGMA_OUT="${PROJECT_DIR}/magma_results/step2_gene_low_mvpa.genes.out"
MAGMA_TISSUE="${PROJECT_DIR}/magma_results/low_mvpa_tissue.gsa.out"
INDEX_FILE="${PROJECT_DIR}/flames_input/indexfile.txt"
FLAMES_OUT_DIR="${PROJECT_DIR}/flames_output"

# ========= VEP & CADD 本地配置 =========
REF_DIR="/home/user2/jzhang_data/Dissertation/code/flames_ref_data"
VEP_CACHE="${REF_DIR}/vep_cache"
CADD_FILE="${REF_DIR}/whole_genome_SNVs.tsv.gz"
TABIX_CMD="${HOME}/miniforge3/envs/FLAMES/bin/tabix"
VEP_CMD="${HOME}/miniforge3/envs/FLAMES/bin/vep"
# ============================

mkdir -p "$FLAMES_OUT_DIR"

# Step 1: Annotate（使用本地 VEP 和 CADD）
echo "[$(date)] Running FLAMES annotate..."
python "${FLAMES_DIR}/FLAMES.py" annotate \
    -a  "$ANNOT_DIR" \
    -p  "$POPS_PREDS" \
    -m  "$MAGMA_OUT" \
    -mt "$MAGMA_TISSUE" \
    -id "$INDEX_FILE" \
    -cv   "$VEP_CMD" \
    -vc  "$VEP_CACHE" \
    -t  "$TABIX_CMD" \
    -cf "$CADD_FILE"

# Step 2: Scoring
echo "[$(date)] Running FLAMES scoring..."
python "${FLAMES_DIR}/FLAMES.py" FLAMES \
    -id "$INDEX_FILE" \
    -o  "$FLAMES_OUT_DIR"

echo "[$(date)] Done. Output: $FLAMES_OUT_DIR"

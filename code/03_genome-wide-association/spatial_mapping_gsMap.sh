#!/bin/bash
set -e
# ============================================================
# 路径与参数
# ============================================================
WORKDIR='/home/user2/jzhang_data/Dissertation/results/high_pa_obesity_uncouple_results/gsmap_results'
ST_DIR='/home/user2/jzhang_data/Dissertation/code/gsMap_example_data/ST/human_embryos'
RESOURCE_DIR='/home/user2/jzhang_data/Dissertation/code/gsMap_resource'
SUMSTATS='/home/Disk/Data2/imputation.bgen/PA_obesity_uncouple_GWAS/high_mvpa_bmi_uncouple_gsMap_format.sumstats.gz'
TRAIT_NAME='high_pa_obesity_uncouple'

# 关键参数（根据实际 h5ad 结构修改）
ANNOTATION='substructure'      # ← obs 中的列名（81 个细分结构）
DATA_LAYER='counts'            # ← layers 中的层名（注意是 counts 不是 count）

SAMPLES=(
    "CS20_E1S6.bin50.substructure"
    "CS20_E1S7.bin50.substructure"
    "CS23_E1S2.bin50.substructure"
    "CS23_E1S3.bin50.substructure"
    "CS23_E2S5.bin50.substructure"
    "CS23_E2S10.bin50.substructure"
)

mkdir -p ${WORKDIR}

# ============================================================
# Step 0: 计算 Slice Mean
# ============================================================
SLICE_MEAN_FILE="${WORKDIR}/sample_slice_mean.parquet"

gsmap create_slice_mean \
    --sample_name_list ${SAMPLES[@]} \
    --h5ad_list \
        "${ST_DIR}/CS20_E1S6.bin50.substructure.h5ad" \
        "${ST_DIR}/CS20_E1S7.bin50.substructure.h5ad" \
        "${ST_DIR}/CS23_E1S2.bin50.substructure.h5ad" \
        "${ST_DIR}/CS23_E1S3.bin50.substructure.h5ad" \
        "${ST_DIR}/CS23_E2S5.bin50.substructure.h5ad" \
        "${ST_DIR}/CS23_E2S10.bin50.substructure.h5ad" \
    --slice_mean_output_file "${SLICE_MEAN_FILE}" \
    --data_layer "${DATA_LAYER}"

# ============================================================
# 对每个样本运行 Step 1~6
# ============================================================
for SAMPLE in "${SAMPLES[@]}"
do
    echo "=========================================="
    echo " Processing sample: ${SAMPLE}"
    echo "=========================================="

    H5AD_PATH="${ST_DIR}/${SAMPLE}.h5ad"

    # ---- Step 1: find_latent_representations ----
    gsmap run_find_latent_representations \
        --workdir "${WORKDIR}" \
        --sample_name "${SAMPLE}" \
        --input_hdf5_path "${H5AD_PATH}" \
        --annotation "${ANNOTATION}" \
        --data_layer "${DATA_LAYER}"

    # ---- Step 2: latent_to_gene （引入 slice mean，人类数据无需 homolog） ----
    gsmap run_latent_to_gene \
        --workdir "${WORKDIR}" \
        --sample_name "${SAMPLE}" \
        --annotation "${ANNOTATION}" \
        --latent_representation 'latent_GVAE' \
        --num_neighbour 51 \
        --num_neighbour_spatial 201 \
        --gM_slices "${SLICE_MEAN_FILE}"

    # ---- Step 3: generate_ldscore ----
    for CHROM in {1..22}
    do
        gsmap run_generate_ldscore \
            --workdir "${WORKDIR}" \
            --sample_name "${SAMPLE}" \
            --chrom ${CHROM} \
            --bfile_root "${RESOURCE_DIR}/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC" \
            --keep_snp_root "${RESOURCE_DIR}/LDSC_resource/hapmap3_snps/hm" \
            --gtf_annotation_file "${RESOURCE_DIR}/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf" \
            --gene_window_size 50000 \
            --enhancer_annotation_file "${RESOURCE_DIR}/genome_annotation/enhancer/by_tissue/ALL/ABC_roadmap_merged.bed" \
            --snp_multiple_enhancer_strategy 'max_mkscore' \
            --gene_window_enhancer_priority 'gene_window_first'
    done

    # ---- Step 4: spatial_ldsc ----
    gsmap run_spatial_ldsc \
        --workdir "${WORKDIR}" \
        --sample_name "${SAMPLE}" \
        --trait_name "${TRAIT_NAME}" \
        --sumstats_file "${SUMSTATS}" \
        --w_file "${RESOURCE_DIR}/LDSC_resource/weights_hm3_no_hla/weights." \
        --num_processes 16

    # ---- Step 5: cauchy_combination （单样本内） ----
    gsmap run_cauchy_combination \
        --workdir "${WORKDIR}" \
        --sample_name "${SAMPLE}" \
        --trait_name "${TRAIT_NAME}" \
        --annotation "${ANNOTATION}"

    # ---- Step 6: 单样本报告 ----
    gsmap run_report \
        --workdir "${WORKDIR}" \
        --sample_name "${SAMPLE}" \
        --trait_name "${TRAIT_NAME}" \
        --annotation "${ANNOTATION}" \
        --sumstats_file "${SUMSTATS}" \
        --top_corr_genes 50
done

# ============================================================
# Step 7: 跨样本聚合 Cauchy combination
# ============================================================
gsmap run_cauchy_combination \
    --workdir "${WORKDIR}" \
    --sample_name_list ${SAMPLES[@]} \
    --trait_name "${TRAIT_NAME}" \
    --annotation "${ANNOTATION}" \
    --output_file "${WORKDIR}/combined_${TRAIT_NAME}_cauchy_combination.csv.gz"




    
# ============================================================
# 路径与参数
# ============================================================
WORKDIR='/home/user2/jzhang_data/Dissertation/results/high_pa_obesity_uncouple_results/gsmap_results'
ST_DIR='/home/user2/jzhang_data/Dissertation/code/gsMap_example_data/ST/human_embryos'
RESOURCE_DIR='/home/user2/jzhang_data/Dissertation/code/gsMap_resource'
SUMSTATS='/home/Disk/Data2/imputation.bgen/PA_obesity_uncouple_GWAS/high_mvpa_bmi_uncouple_gsMap_format.sumstats.gz'
TRAIT_NAME='high_pa_obesity_uncouple'

ANNOTATION='substructure'
DATA_LAYER='counts'

SAMPLE="CS23_E2S10.bin50.substructure"
H5AD_PATH="${ST_DIR}/${SAMPLE}.h5ad"
SLICE_MEAN_FILE="${WORKDIR}/sample_slice_mean.parquet"

mkdir -p ${WORKDIR}

# ============================================================
# Step 0: 跳过（已有 slice mean）
# ============================================================
echo "Slice mean already exists: ${SLICE_MEAN_FILE}"

# ============================================================
# Step 1: find_latent_representations
# ============================================================
echo "=========================================="
echo " Processing sample: ${SAMPLE}"
echo "=========================================="

gsmap run_find_latent_representations \
    --workdir "${WORKDIR}" \
    --sample_name "${SAMPLE}" \
    --input_hdf5_path "${H5AD_PATH}" \
    --annotation "${ANNOTATION}" \
    --data_layer "${DATA_LAYER}"

# ============================================================
# Step 2: latent_to_gene
# ============================================================
gsmap run_latent_to_gene \
    --workdir "${WORKDIR}" \
    --sample_name "${SAMPLE}" \
    --annotation "${ANNOTATION}" \
    --latent_representation 'latent_GVAE' \
    --num_neighbour 51 \
    --num_neighbour_spatial 201 \
    --gM_slices "${SLICE_MEAN_FILE}"

# ============================================================
# Step 3: generate_ldscore
# ============================================================
for CHROM in {1..22}
do
    gsmap run_generate_ldscore \
        --workdir "${WORKDIR}" \
        --sample_name "${SAMPLE}" \
        --chrom ${CHROM} \
        --bfile_root "${RESOURCE_DIR}/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC" \
        --keep_snp_root "${RESOURCE_DIR}/LDSC_resource/hapmap3_snps/hm" \
        --gtf_annotation_file "${RESOURCE_DIR}/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf" \
        --gene_window_size 50000 \
        --enhancer_annotation_file "${RESOURCE_DIR}/genome_annotation/enhancer/by_tissue/ALL/ABC_roadmap_merged.bed" \
        --snp_multiple_enhancer_strategy 'max_mkscore' \
        --gene_window_enhancer_priority 'gene_window_first'
done

# ============================================================
# Step 4: spatial_ldsc
# ============================================================
gsmap run_spatial_ldsc \
    --workdir "${WORKDIR}" \
    --sample_name "${SAMPLE}" \
    --trait_name "${TRAIT_NAME}" \
    --sumstats_file "${SUMSTATS}" \
    --w_file "${RESOURCE_DIR}/LDSC_resource/weights_hm3_no_hla/weights." \
    --num_processes 32

# ============================================================
# Step 5: cauchy_combination
# ============================================================
gsmap run_cauchy_combination \
    --workdir "${WORKDIR}" \
    --sample_name "${SAMPLE}" \
    --trait_name "${TRAIT_NAME}" \
    --annotation "${ANNOTATION}"

# ============================================================
# Step 6: 报告
# ============================================================
gsmap run_report \
    --workdir "${WORKDIR}" \
    --sample_name "${SAMPLE}" \
    --trait_name "${TRAIT_NAME}" \
    --annotation "${ANNOTATION}" \
    --sumstats_file "${SUMSTATS}" \
    --top_corr_genes 50

echo "=========================================="
echo " Done: ${SAMPLE}"
echo "=========================================="

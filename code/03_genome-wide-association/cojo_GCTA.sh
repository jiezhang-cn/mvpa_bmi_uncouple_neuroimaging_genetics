#GCTA COJO Conditional Analysis Pipeline
# ─────────────────────────────────────────────────────────────────────────────
# Step 1  : 逐染色体 GCTA --cojo-slct (直接用 per-chr LD reference)
# Step 2  : 合并所有染色体结果
# Step 2.5: 折叠 long-range high-LD 区域内的多余信号
# Step 3  : 过滤 (原始 GWAS P < 5e-8)
#
# Reference for long-range LD regions:
#   Price AL et al. Am J Hum Genet. 2008;83(1):132-135.
#   Anderson CA et al. Nat Protoc. 2010;5(9):1564-73.
###############################################################################

set -e

# ============================================================================
# 配置参数
# ============================================================================

GCTA="$HOME/miniforge3/envs/gcta_env/bin/gcta"


# LD 参考面板目录 (per-chromosome: ref_chr1, ref_chr2, ..., ref_chr22)
REF_DIR="/test/jzhang_data_depository/gcta_ukb_unrelated_ref"

# GWAS summary statistics
GWAS_DIR="/home/Disk/Data2/imputation.bgen/PA_obesity_uncouple_GWAS"

# GWAS 文件 -> 对应输出目录
declare -A GWAS_MAP
GWAS_MAP["${GWAS_DIR}/high_mvpa_bmi_uncouple.ma"]="/home/user2/jzhang_data/Dissertation/results/high_pa_obesity_uncouple_results/gcta_cojo_results"
GWAS_MAP["${GWAS_DIR}/low_mvpa_bmi_uncouple.ma"]="/home/user2/jzhang_data/Dissertation/results/low_pa_obesity_uncouple_results/gcta_cojo_results"

# 按顺序执行的文件列表
GWAS_FILES=(
    "${GWAS_DIR}/high_mvpa_bmi_uncouple.ma"
    "${GWAS_DIR}/low_mvpa_bmi_uncouple.ma"
)

# 分析参数
COJO_P="5e-8"
ORIG_P="5e-8"
MAF=0.01
COJO_WIND=10000
COJO_COLLINEAR=0.9
THREADS=24

# ============================================================================
# Long-range high-LD 区域定义 (hg19 / GRCh37)
# ============================================================================
# 格式: "CHR START_BP END_BP REGION_NAME"
# 来源:
#   - Price AL et al. Am J Hum Genet 2008;83:132-135
#   - Anderson CA et al. Nat Protoc 2010;5:1564-73
#   - 部分区域边界略有扩展以覆盖不同研究中报告的范围

LONG_RANGE_LD_REGIONS=(
    "1   48000000   52000000   1p33.3-p32.3"
    "2   86000000  100500000   2p11-q12.1_LCT"
    "2  134500000  138000000   2q21.1"
    "2  183000000  190000000   2q32.2-q33.1"
    "3   47500000   50000000   3p21.1"
    "3   83500000   87000000   3p12.1"
    "3   89000000   97500000   3p12.1-p11.1"
    "5   44000000   51500000   5p12-p11"
    "5   98000000  100500000   5q21.1"
    "5  129000000  132000000   5q23.3-q31.1"
    "5  135500000  138500000   5q31.1-q31.2"
    "6   25500000   33500000   6p22.1-p21.1_MHC"
    "6   57000000   64000000   6p12.3-q12"
    "6  140000000  142500000   6q24.1"
    "7   55000000   66000000   7p11.2-q11.23"
    "8    8000000   12000000   8p23.1_inversion"
    "8   43000000   50000000   8p12-p11.1"
    "8  112000000  115000000   8q23.3-q24.11"
    "10  37000000   43000000   10p11.21-q11.1"
    "11  46000000   57000000   11p11.2-q13.1"
    "11  87500000   90500000   11q14.1-q14.3"
    "12  33000000   40000000   12p11.22-q12"
    "12 109500000  112000000   12q24.12-q24.13"
    "20  32000000   34500000   20q11.21-q11.23"
)

# ============================================================================
# 函数
# ============================================================================

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# ============================================================================
# Step 1: GCTA COJO per chromosome (直接用 per-chr bfile)
# ============================================================================

run_cojo_per_chr() {
    local GWAS_FILE=$1
    local TRAIT=$2
    local OUT_DIR=$3

    for CHR in $(seq 1 22); do

        local BFILE="${REF_DIR}/ref_chr${CHR}"
        local OUT_PREFIX="${OUT_DIR}/${TRAIT}_chr${CHR}"

        # 检查该染色体参考文件是否存在
        if [ ! -f "${BFILE}.bed" ] || [ ! -f "${BFILE}.bim" ] || [ ! -f "${BFILE}.fam" ]; then
            log_msg "    [WARNING] Chr${CHR} reference files incomplete, skipping."
            continue
        fi

        log_msg "  Running COJO on Chr${CHR}..."

        ${GCTA} \
            --bfile "${BFILE}" \
            --maf ${MAF} \
            --cojo-file "${GWAS_FILE}" \
            --cojo-slct \
            --cojo-p ${COJO_P} \
            --cojo-wind ${COJO_WIND} \
            --cojo-collinear ${COJO_COLLINEAR} \
            --thread-num ${THREADS} \
            --out "${OUT_PREFIX}" \
            > "${OUT_PREFIX}.log" 2>&1

        if [ -f "${OUT_PREFIX}.jma.cojo" ]; then
            local N_SIG=$(tail -n +2 "${OUT_PREFIX}.jma.cojo" | wc -l)
            log_msg "    Chr${CHR}: ${N_SIG} independent signal(s)"
        else
            log_msg "    Chr${CHR}: No independent signals"
        fi

    done
}

# ============================================================================
# Step 2: Merge results across chromosomes
# ============================================================================

merge_results() {
    local TRAIT=$1
    local OUT_DIR=$2
    local MERGED_FILE="${OUT_DIR}/${TRAIT}_all_chr.jma.cojo"

    local HEADER_DONE=false

    for CHR in $(seq 1 22); do
        local RESULT="${OUT_DIR}/${TRAIT}_chr${CHR}.jma.cojo"
        if [ -f "${RESULT}" ]; then
            if [ "${HEADER_DONE}" = false ]; then
                cat "${RESULT}" > "${MERGED_FILE}"
                HEADER_DONE=true
            else
                tail -n +2 "${RESULT}" >> "${MERGED_FILE}"
            fi
        fi
    done

    if [ "${HEADER_DONE}" = true ]; then
        local N_TOTAL=$(tail -n +2 "${MERGED_FILE}" | wc -l)
        log_msg "  Merged results: ${MERGED_FILE}"
        log_msg "  Total conditional independent signals: ${N_TOTAL}"
        echo "${N_TOTAL}"
    else
        log_msg "  No significant results across all chromosomes."
        echo "0"
    fi
}

# ============================================================================
# Step 2.5: 折叠 long-range high-LD 区域内的多余信号
# ============================================================================

collapse_long_range_ld() {
    local TRAIT=$1
    local OUT_DIR=$2
    local INPUT_FILE="${OUT_DIR}/${TRAIT}_all_chr.jma.cojo"
    local OUTPUT_FILE="${OUT_DIR}/${TRAIT}_all_chr_lrcollapsed.jma.cojo"
    local REMOVED_FILE="${OUT_DIR}/${TRAIT}_all_chr_lr_removed.txt"

    if [ ! -f "${INPUT_FILE}" ]; then
        log_msg "  [WARNING] Input file not found for LR-LD collapsing."
        return
    fi

    local N_BEFORE=$(tail -n +2 "${INPUT_FILE}" | wc -l)

    # ── 生成区域定义的临时文件 ──
    local REGION_FILE="${OUT_DIR}/.tmp_lr_ld_regions.txt"
    > "${REGION_FILE}"
    for REGION in "${LONG_RANGE_LD_REGIONS[@]}"; do
        echo "${REGION}" >> "${REGION_FILE}"
    done

    # ── AWK 实现折叠逻辑 ──
    awk -v region_file="${REGION_FILE}" -v removed_file="${REMOVED_FILE}" '
    BEGIN {
        n_regions = 0
        while ((getline line < region_file) > 0) {
            n_regions++
            split(line, f)
            reg_chr[n_regions]   = f[1] + 0
            reg_start[n_regions] = f[2] + 0
            reg_end[n_regions]   = f[3] + 0
            reg_name[n_regions]  = f[4]
        }
        close(region_file)
    }

    NR == 1 {
        header = $0
        next
    }
    {
        chr = $1 + 0
        snp = $2
        bp  = $3 + 0
        pJ  = $13 + 0

        assigned_region = 0
        for (i = 1; i <= n_regions; i++) {
            if (chr == reg_chr[i] && bp >= reg_start[i] && bp <= reg_end[i]) {
                assigned_region = i
                break
            }
        }

        if (assigned_region == 0) {
            keep[NR]    = 1
            lines[NR]   = $0
            region[NR]  = 0
        } else {
            keep[NR]    = 1
            lines[NR]   = $0
            region[NR]  = assigned_region
            snps[NR]    = snp
            pJs[NR]     = pJ

            if (!(assigned_region in best_pJ) || pJ < best_pJ[assigned_region]) {
                best_pJ[assigned_region]  = pJ
                best_idx[assigned_region] = NR
            }

            region_members[assigned_region] = region_members[assigned_region] " " NR
        }
    }

    END {
        for (r = 1; r <= n_regions; r++) {
            if (!(r in region_members)) continue

            n_mem = split(region_members[r], members)
            if (n_mem <= 1) continue

            for (j = 1; j <= n_mem; j++) {
                idx = members[j] + 0
                if (idx == 0) continue
                if (idx != best_idx[r]) {
                    keep[idx] = 0
                }
            }
        }

        print header
        for (i = 2; i <= NR; i++) {
            if (keep[i] == 1 && lines[i] != "") {
                print lines[i]
            }
        }

        printf "SNP\tChr\tbp\tpJ\tRegion\tKept_SNP\n" > removed_file
        for (r = 1; r <= n_regions; r++) {
            if (!(r in region_members)) continue
            n_mem = split(region_members[r], members)
            if (n_mem <= 1) continue

            for (j = 1; j <= n_mem; j++) {
                idx = members[j] + 0
                if (idx == 0) continue
                if (keep[idx] == 0) {
                    split(lines[idx], f)
                    split(lines[best_idx[r]], bf)
                    printf "%s\t%s\t%s\t%s\t%s\t%s\n", \
                        f[2], f[1], f[3], f[13], reg_name[r], bf[2] > removed_file
                }
            }
        }
        close(removed_file)
    }
    ' "${INPUT_FILE}" > "${OUTPUT_FILE}"

    local N_AFTER=$(tail -n +2 "${OUTPUT_FILE}" | wc -l)
    local N_REMOVED=$((N_BEFORE - N_AFTER))

    log_msg "  Long-range LD collapsing (hg19 coordinates):"
    log_msg "    Before:  ${N_BEFORE} signals"
    log_msg "    Removed: ${N_REMOVED} signals (redundant within LR-LD regions)"
    log_msg "    After:   ${N_AFTER} signals"

    if [ ${N_REMOVED} -gt 0 ] && [ -f "${REMOVED_FILE}" ]; then
        log_msg "    Removed signals detail: ${REMOVED_FILE}"
        log_msg "    ──────────────────────────────────────"
        while IFS=$'\t' read -r SNP CHR BP PJ REGION KEPT; do
            if [ "${SNP}" = "SNP" ]; then continue; fi
            log_msg "    Removed: ${SNP} (chr${CHR}:${BP}, pJ=${PJ})"
            log_msg "             Region: ${REGION}, kept: ${KEPT}"
        done < "${REMOVED_FILE}"
        log_msg "    ──────────────────────────────────────"
    fi

    rm -f "${REGION_FILE}"

    echo "${N_AFTER}"
}

# ============================================================================
# Step 3: 过滤 (原始 GWAS P < threshold)
# ============================================================================

filter_by_original_p() {
    local GWAS_FILE=$1
    local TRAIT=$2
    local OUT_DIR=$3

    local INPUT_FILE="${OUT_DIR}/${TRAIT}_all_chr_lrcollapsed.jma.cojo"
    if [ ! -f "${INPUT_FILE}" ]; then
        INPUT_FILE="${OUT_DIR}/${TRAIT}_all_chr.jma.cojo"
    fi

    local FILTERED_FILE="${OUT_DIR}/${TRAIT}_all_chr_filtered.jma.cojo"

    if [ ! -f "${INPUT_FILE}" ]; then
        log_msg "  [WARNING] Input file not found, skipping filter step."
        return
    fi

    awk -v pthresh="${ORIG_P}" '
    BEGIN { pthresh_num = pthresh + 0 }
    FNR == NR {
        if (FNR == 1) next
        orig_p[$1] = $7 + 0
        next
    }
    FNR != NR {
        if (FNR == 1) { print $0; next }
        snp = $2
        if (snp in orig_p && orig_p[snp] < pthresh_num) {
            print $0
        }
    }
    ' "${GWAS_FILE}" "${INPUT_FILE}" > "${FILTERED_FILE}"

    local N_FILTERED=$(tail -n +2 "${FILTERED_FILE}" | wc -l)
    log_msg "  After filtering (original P < ${ORIG_P}): ${N_FILTERED} index SNPs"
    log_msg "  Final output: ${FILTERED_FILE}"
}

# ============================================================================
# 主程序
# ============================================================================

main() {
    log_msg "============================================================"
    log_msg " GCTA COJO Conditional Analysis Pipeline"
    log_msg " GCTA:                    ${GCTA}"
    log_msg " LD reference (per-chr):  ${REF_DIR}/ref_chr{1..22}"
    log_msg " COJO P threshold:        ${COJO_P}"
    log_msg " Original GWAS P filter:  ${ORIG_P}"
    log_msg " Genome build:            hg19 (GRCh37)"
    log_msg " Long-range LD regions:   ${#LONG_RANGE_LD_REGIONS[@]}"
    log_msg "============================================================"
    echo ""

    # 检查 GCTA
    if [ ! -x "${GCTA}" ]; then
        log_msg "[ERROR] GCTA not found or not executable: ${GCTA}"
        exit 1
    fi
    log_msg "GCTA: $(${GCTA} --version 2>&1 | head -1 || echo 'version check failed')"

    # 检查参考面板
    local REF_COUNT=0
    for CHR in $(seq 1 22); do
        if [ -f "${REF_DIR}/ref_chr${CHR}.bed" ]; then
            REF_COUNT=$((REF_COUNT + 1))
        fi
    done
    log_msg "LD reference: ${REF_COUNT}/22 chromosomes found"

    if [ ${REF_COUNT} -eq 0 ]; then
        log_msg "[ERROR] No reference panel files found in ${REF_DIR}!"
        exit 1
    fi

    # 打印 long-range LD 区域摘要
    log_msg ""
    log_msg "Long-range high-LD regions (hg19):"
    for REGION in "${LONG_RANGE_LD_REGIONS[@]}"; do
        read -r CHR START END NAME <<< "${REGION}"
        local SIZE_MB=$(echo "scale=1; (${END} - ${START}) / 1000000" | bc)
        log_msg "  chr${CHR}:${START}-${END} (${SIZE_MB} Mb) ${NAME}"
    done
    echo ""

    # ── 逐个 GWAS 文件处理 ──
    for GWAS_FILE in "${GWAS_FILES[@]}"; do

        TRAIT=$(basename "${GWAS_FILE}" .ma)
        OUT_DIR="${GWAS_MAP[${GWAS_FILE}]}"

        mkdir -p "${OUT_DIR}"

        echo ""
        log_msg "========================================================"
        log_msg " Trait:  ${TRAIT}"
        log_msg " Input:  ${GWAS_FILE}"
        log_msg " Output: ${OUT_DIR}"
        log_msg "========================================================"

        if [ ! -f "${GWAS_FILE}" ]; then
            log_msg "[ERROR] GWAS file not found: ${GWAS_FILE}"
            continue
        fi

        log_msg "  .ma header: $(head -1 "${GWAS_FILE}")"
        log_msg "  .ma lines:  $(wc -l < "${GWAS_FILE}")"

        # Step 1
        log_msg ""
        log_msg "── Step 1: GCTA --cojo-slct per chromosome ──"
        run_cojo_per_chr "${GWAS_FILE}" "${TRAIT}" "${OUT_DIR}"

        # Step 2
        log_msg ""
        log_msg "── Step 2: Merge results across chromosomes ──"
        N_TOTAL=$(merge_results "${TRAIT}" "${OUT_DIR}")

        if [ "${N_TOTAL}" = "0" ]; then
            log_msg "  No signals to process, moving to next trait."
            continue
        fi

        # Step 2.5
        log_msg ""
        log_msg "── Step 2.5: Collapse signals in long-range high-LD regions ──"
        log_msg "  Ref: Price AL et al. Am J Hum Genet 2008;83:132-135"
        N_AFTER_LR=$(collapse_long_range_ld "${TRAIT}" "${OUT_DIR}")

        if [ "${N_AFTER_LR}" = "0" ]; then
            log_msg "  No signals remaining after LR-LD collapsing."
            continue
        fi

        # Step 3
        log_msg ""
        log_msg "── Step 3: Filter by original GWAS P < ${ORIG_P} ──"
        filter_by_original_p "${GWAS_FILE}" "${TRAIT}" "${OUT_DIR}"

    done

    echo ""
    log_msg "============================================================"
    log_msg " All analyses completed!"
    log_msg "============================================================"
}

# 日志目录
LOG_DIR="/home/user2/jzhang_data/Dissertation/code"
mkdir -p "${LOG_DIR}"

# 运行并保存日志
main 2>&1 | tee "${LOG_DIR}/gcta_cojo_pipeline_$(date '+%Y%m%d_%H%M%S').log"

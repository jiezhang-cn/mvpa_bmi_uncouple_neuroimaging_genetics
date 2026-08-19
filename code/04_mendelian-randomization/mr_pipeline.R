#!/usr/bin/env Rscript
###############################################################################
#  Bi-directional MR + CAUSE
#  Forward : Psychiatric disorders -> MVPA_{high,low}
#  Reverse : MVPA_{high,low}        -> Psychiatric disorders
#
#  Changes (this revision):
#   1. Resume: skip pairs whose output files already exist (and are non-empty)
#   2. Timeout protection on CAUSE (no thread limiting); uses ELAPSED time only
#   3. STRICT pval_thresh (1e-4) for PTSD / AUDIT_P / ANX
#   4. Fresh run.log per invocation (overwrite, not append)
#   5. 2SMR auto-retry with p<1e-5 IVs if <3 SNPs after harmonise
#   6. CAUSE drops extreme SNPs (|z|>38, SE bad, non-finite) post-merge
#   7. Stronger resume — pre-check both directions BEFORE any heavy IO
#   8. ★ NEW: posterior SNP hard cap (POSTERIOR_SNP_CAP) to avoid cause() blow-up
#   9. ★ NEW: withTimeout uses cpu=Inf, elapsed=... (wall-clock only)
###############################################################################
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr)
  library(TwoSampleMR); library(MRPRESSO); library(ieugwasr)
  library(cause); library(ggplot2); library(R.utils)
})

# ---------------- Parameters ------------------------------------------------
op <- list(
  make_option("--out"),       make_option("--mvpa_high"), make_option("--mvpa_low"),
  make_option("--gwas_dir"),  make_option("--cd_dir"),    make_option("--prep_dir"),
  make_option("--bfile"),     make_option("--plink_bin")
)
opt <- parse_args(OptionParser(option_list=op))

dir.create(file.path(opt$out,"cache"),         showWarnings=FALSE)
dir.create(file.path(opt$out,"plots","cause"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out,"results"),       showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out,"sensitivity"),   showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out,"harmonised"),    showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out,"exposure_dat"),  showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out,"plots"),         showWarnings=FALSE, recursive=TRUE)

# ============================================================================
# 0. Logging  (message sink cannot use split=TRUE)
# ============================================================================
log_path <- file.path(opt$out, "run.log")
log_con  <- file(log_path, open = "wt")   # overwrite

sink(log_con, split = TRUE,  type = "output")
sink(log_con, split = FALSE, type = "message")

cat("=== Run started:", format(Sys.time()), "===\n")
cat("Log file:", log_path, "\n\n")

on.exit({
  cat("\n=== Run finished:", format(Sys.time()), "===\n")
  try(sink(type = "message"), silent = TRUE)
  try(sink(),                 silent = TRUE)
  try(close(log_con),         silent = TRUE)
}, add = TRUE)

PLINK <- opt$plink_bin
BFILE <- opt$bfile

# ---- timing limits (wall-clock, seconds) ----
CAUSE_TIMEOUT    <- 7200   # 120 min for cause()
NUISANCE_TIMEOUT <- 5400   # 90 min for est_cause_params()

# ---- ★ disorders requiring stricter exposure-pval threshold in CAUSE ----
STRICT_DISORDERS <- c("PTSD", "AUDIT_P", "ANX")

# ---- ★ posterior SNP hard cap (CAUSE recommends ~1000) ----
POSTERIOR_SNP_CAP <- 2000

# Helper: file exists AND has non-zero size
file_done <- function(p) {
  isTRUE(file.exists(p) && file.info(p)$size > 0)
}

# ============================================================================
# 1. GWAS column mappings  (OUD removed per user request)
# ============================================================================
mvpa_list <- list(
  MVPA_high = list(file=opt$mvpa_high,
                   snp="SNP", chr="CHR", bp="BP",
                   ea="ALLELE1", oa="ALLELE0", eaf="A1FREQ", info="INFO",
                   beta="BETA", se="SE", pval="P_BOLT_LMM",
                   n=106037),
  MVPA_low  = list(file=opt$mvpa_low,
                   snp="SNP", chr="CHR", bp="BP",
                   ea="ALLELE1", oa="ALLELE0", eaf="A1FREQ", info="INFO",
                   beta="BETA", se="SE", pval="P_BOLT_LMM",
                   n=105222)
)

disorders <- list(
  ## ---- 6 cross-disorder factors (Grotzinger 2025) ----
  F1_Compulsive    = list(file=file.path(opt$cd_dir,"F1_CompulsiveDisorders_2025.tsv.gz"),
                          snp="SNP",beta="BETA",se="SE",ea="A1",oa="A2",pval="P",n=54100),
  F2_SCZBIP        = list(file=file.path(opt$cd_dir,"F2_SchizophreniaBipolar_2025.tsv.gz"),
                          snp="SNP",beta="BETA",se="SE",ea="A1",oa="A2",pval="P",n=127202),
  F3_Neurodev      = list(file=file.path(opt$cd_dir,"F3_Neurodevelopmental_2025.tsv.gz"),
                          snp="SNP",beta="BETA",se="SE",ea="A1",oa="A2",pval="P",n=84760),
  F4_Internalizing = list(file=file.path(opt$cd_dir,"F4_Internalizing_2025.tsv.gz"),
                          snp="SNP",beta="BETA",se="SE",ea="A1",oa="A2",pval="P",n=1637337),
  F5_SubstanceUse  = list(file=file.path(opt$cd_dir,"F5_SubstanceUse_2025.tsv.gz"),
                          snp="SNP",beta="BETA",se="SE",ea="A1",oa="A2",pval="P",n=313395),
  PFactor          = list(file=file.path(opt$cd_dir,"PFactor_2025.tsv.gz"),
                          snp="SNP",beta="BETA",se="SE",ea="A1",oa="A2",pval="P",n=2168621),

  ## ---- Individual disorders ----
  AN    = list(file=file.path(opt$prep_dir,"AN.tsv.gz"),
               snp="ID",beta="BETA",se="SE",ea="ALT",oa="REF",pval="PVAL",
               info="IMPINFO", ncase=16992,ncontrol=55525),

  OCD   = list(file=file.path(opt$gwas_dir,"daner_OCDmeta_wo23andMe_LOOUKBB_080425.gz"),
               snp="SNP",or="OR",se="SE",ea="A1",oa="A2",
               eaf="FRQ_A_22717", info="INFO", pval="P",
               ncase=22717, ncontrol=988884),

  ANX   = list(file=file.path(opt$gwas_dir,"ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz"),
               snp="SNP",or="OR",se="SE",ea="A1",oa="A2",
               eaf="FRQ_A_122083", info="INFO", pval="P",
               ncase=122083, ncontrol=729602),

  SCZ   = list(file=file.path(opt$prep_dir,"SCZ.tsv.gz"),
               snp="ID",beta="BETA",se="SE",ea="A1",oa="A2",
               eaf="FCAS", info="IMPINFO", pval="PVAL", n_col="NEFF"),

  BIP   = list(file=file.path(opt$gwas_dir,"daner_bip_pgc3_nm_noukbiobank.gz"),
               snp="SNP",or="OR",se="SE",ea="A1",oa="A2",
               eaf="FRQ_A_40463", info="INFO", pval="P",
               ncase=40463, ncontrol=313436),

  ASD   = list(file=file.path(opt$gwas_dir,"iPSYCH-PGC_ASD_Nov2017.gz"),
               snp="SNP",or="OR",se="SE",ea="A1",oa="A2",
               info="INFO", pval="P", n=46350),

  ADHD  = list(file=file.path(opt$gwas_dir,"ADHD2022_iPSYCH_deCODE_PGC.meta.gz"),
               snp="SNP",or="OR",se="SE",ea="A1",oa="A2",
               eaf="FRQ_A_38691", info="INFO", pval="P",
               ncase=38691, ncontrol=186843),

  PTSD  = list(file=file.path(opt$prep_dir,"PTSD.tsv.gz"),
               snp="ID",z="Z",ea="A1",oa="A2",
               eaf="FREQ", pval="P", n_col="NEFF"),

  MDD   = list(file=file.path(opt$gwas_dir,"daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz"),
               snp="SNP",or="OR",se="SE",ea="A1",oa="A2",
               eaf="FRQ_A_45396", info="INFO", pval="P",
               ncase=45396, ncontrol=97250),

  CUD   = list(file=file.path(opt$gwas_dir,"CUD_EUR_full_public_11.14.2020.gz"),
               snp="SNP",z="Z",ea="A1",oa="A2",pval="P",n_col="N"),

  AUDIT_P = list(file=file.path(opt$prep_dir,"AUDIT_P.tsv.gz"),
                 snp="SNP",beta="BETA",se="SE",ea="A1",oa="A2",
                 info="INFO", pval="P", n_col="N")
)

# ============================================================================
# 2. read_gwas - safe version
# ============================================================================
read_gwas <- function(meta, name, cache_dir){
  fcache <- file.path(cache_dir, paste0(name,".rds"))
  if (file.exists(fcache)) {
    cat("  [cache]", name, "\n")
    return(readRDS(fcache))
  }
  cat("  reading", name, "...\n")
  dt <- fread(meta$file, data.table=FALSE)

  pick <- function(col, type=c("num","chr")){
    type <- match.arg(type)
    if (is.null(col)) return(rep(NA, nrow(dt)))
    if (!col %in% names(dt)){
      cat(sprintf("    [warn] column '%s' NOT in %s -> filled NA\n", col, name))
      return(rep(NA, nrow(dt)))
    }
    v <- dt[[col]]
    if (type=="num") suppressWarnings(as.numeric(v)) else as.character(v)
  }

  out <- data.frame(
    SNP            = pick(meta$snp,  "chr"),
    effect_allele  = toupper(pick(meta$ea, "chr")),
    other_allele   = toupper(pick(meta$oa, "chr")),
    pval           = pick(meta$pval, "num"),
    stringsAsFactors = FALSE
  )

  if (!is.null(meta$beta)) {
    out$beta <- pick(meta$beta, "num")
    out$se   <- pick(meta$se,   "num")
  } else if (!is.null(meta$or)) {
    or_v <- pick(meta$or, "num")
    out$beta <- log(or_v)
    out$se   <- pick(meta$se, "num")
  } else if (!is.null(meta$z)) {
    out$beta <- pick(meta$z, "num")
    out$se   <- rep(1, nrow(out))
  }

  out$eaf  <- pick(meta$eaf,  "num")
  out$info <- pick(meta$info, "num")

  if (!is.null(meta$n) && is.numeric(meta$n) && length(meta$n)==1) {
    out$samplesize <- meta$n
  } else if (!is.null(meta$n_col)) {
    out$samplesize <- pick(meta$n_col, "num")
  } else if (!is.null(meta$ncase) && !is.null(meta$ncontrol)) {
    out$ncase      <- meta$ncase
    out$ncontrol   <- meta$ncontrol
    out$samplesize <- meta$ncase + meta$ncontrol
  } else if (!is.null(meta$ncase_col)) {
    out$ncase      <- pick(meta$ncase_col,    "num")
    out$ncontrol   <- pick(meta$ncontrol_col, "num")
    out$samplesize <- out$ncase + out$ncontrol
  } else {
    out$samplesize <- NA_real_
  }

  keep <- !is.na(out$beta) & !is.na(out$se) & !is.na(out$pval) &
          out$se > 0 & is.finite(out$beta) &
          !is.na(out$SNP) & !is.na(out$effect_allele) & !is.na(out$other_allele)
  out <- out[keep, , drop=FALSE]

  if (!is.null(meta$info) && meta$info %in% names(dt)) {
    out <- out[is.na(out$info) | out$info >= 0.8, , drop=FALSE]
  }

  out$Phenotype <- name
  cat(sprintf("    -> %d SNPs kept\n", nrow(out)))
  saveRDS(out, fcache)
  out
}

# ============================================================================
# 3. IV selection
# ============================================================================
get_iv <- function(dat, p=5e-8, r2=0.001, kb=10000){
  sig <- dat[!is.na(dat$pval) & dat$pval < p, ]
  if (nrow(sig) < 3) {
    cat("    [warn] <3 SNPs at p<", p, "\n")
    return(NULL)
  }
  clp <- tryCatch(
    ld_clump(dplyr::tibble(rsid=sig$SNP, pval=sig$pval),
             clump_r2=r2, clump_kb=kb, plink_bin=PLINK, bfile=BFILE),
    error=function(e){ cat("    clump error:", conditionMessage(e), "\n"); NULL })
  if (is.null(clp) || nrow(clp)==0) return(NULL)
  sig[sig$SNP %in% clp$rsid, ]
}

# Wrapper: try 5e-8 first, then fallback to 1e-5
get_iv_with_fallback <- function(dat, name, cache_path = NULL){
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cat("    [cache IV]", name, "\n")
    x <- tryCatch(fread(cache_path, data.table=FALSE), error = function(e) NULL)
    if (!is.null(x) && nrow(x) >= 3) return(x)
    cat("    cached IV insufficient -> re-extract\n")
  }
  x <- get_iv(dat, p = 5e-8)
  if (is.null(x) || nrow(x) < 3) {
    cat("    relax threshold to p<1e-5\n")
    x <- get_iv(dat, p = 1e-5)
  }
  if (!is.null(x) && !is.null(cache_path)) fwrite(x, cache_path, sep="\t")
  x
}

# ============================================================================
# 4. F-statistic
# ============================================================================
f_stat <- function(d){
  beta <- d$beta.exposure
  se   <- d$se.exposure
  eaf  <- d$eaf.exposure
  n    <- d$samplesize.exposure
  k    <- nrow(d)

  ok <- !is.na(eaf) & !is.na(beta) & !is.na(se) & !is.na(n) &
        eaf > 0 & eaf < 1 & se > 0 & n > 0
  if (sum(ok) > 0) {
    eafx <- eaf[ok]; betax <- beta[ok]; sex <- se[ok]; nx <- n[ok]
    r2_snp <- (2 * eafx * (1-eafx) * betax^2) /
              (2 * eafx * (1-eafx) * betax^2 +
               2 * nx * eafx * (1-eafx) * sex^2)
    R2     <- min(sum(r2_snp, na.rm=TRUE), 0.999)
    n_mean <- mean(nx, na.rm=TRUE)
    Fstat  <- (n_mean - k - 1) / k * R2 / (1 - R2)
    return(data.frame(R2=R2, F=Fstat, k=k))
  }

  ok2 <- is.finite(beta) & is.finite(se) & se > 0
  if (sum(ok2) == 0) return(data.frame(R2=NA_real_, F=NA_real_, k=k))
  betax <- beta[ok2]; sex <- se[ok2]; nx <- n[ok2]
  Fi <- (betax / sex)^2

  if (any(!is.na(nx) & nx > 0)) {
    n_mean <- mean(nx[!is.na(nx) & nx > 0])
    nx[is.na(nx) | nx <= 0] <- n_mean
    r2_snp <- Fi / (Fi + nx - 2)
    R2 <- min(sum(r2_snp, na.rm=TRUE), 0.999)
    Fstat <- (n_mean - k - 1) / k * R2 / (1 - R2)
    return(data.frame(R2=R2, F=Fstat, k=k))
  }

  data.frame(R2=NA_real_, F=mean(Fi), k=k)
}

fmt_num <- function(x, d) ifelse(is.na(x), "NA", formatC(x, digits=d, format="f"))

# ============================================================================
# 4b. Clean extreme SNPs from CAUSE merged dataset
# ============================================================================
clean_cause_X <- function(X){
  if (is.null(X) || nrow(X) == 0) return(X)
  z1 <- X$beta_hat_1 / X$seb1
  z2 <- X$beta_hat_2 / X$seb2
  keep <- is.finite(X$beta_hat_1) & is.finite(X$beta_hat_2) &
          is.finite(X$seb1) & is.finite(X$seb2) &
          X$seb1 > 1e-4 & X$seb2 > 1e-4 &
          X$seb1 < 10   & X$seb2 < 10   &
          is.finite(z1) & is.finite(z2) &
          abs(z1) < 38  & abs(z2) < 38
  n_drop <- sum(!keep, na.rm = TRUE)
  if (n_drop > 0)
    cat(sprintf("    cleaned %d extreme SNPs (|z|>38 / SE bad / non-finite)\n", n_drop))
  X[keep, , drop = FALSE]
}

# ============================================================================
# 5. CAUSE
#    Resume: skip if results/{tag}_cause.rds already exists (non-empty)
#    Timeout: wrap cause() and est_cause_params() in withTimeout()
#             - uses cpu=Inf, elapsed=T (wall-clock only, avoids multithread CPU blowup)
#    STRICT pval_thresh = 1e-4 for PTSD / AUDIT_P / ANX (others 1e-3)
#    Posterior SNP cap: trim to POSTERIOR_SNP_CAP after LD pruning if needed
# ============================================================================
run_cause <- function(exp_full, out_full, exp_name, out_name, outdir){
  tag <- paste0(exp_name, "_to_", out_name)
  cat("  [CAUSE]", tag, "\n")

  # ---- RESUME ----
  rds_path <- file.path(outdir, "results", paste0(tag,"_cause.rds"))
  if (file_done(rds_path)) {
    cat("    [skip] already done ->", rds_path, "\n")
    res <- tryCatch(readRDS(rds_path), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    s <- tryCatch(summary(res, ci_size = 0.95), error = function(e) NULL)
    if (is.null(s)) return(NULL)
    return(data.frame(
      exposure   = exp_name, outcome = out_name,
      n_variants = NA_integer_,
      gamma_est  = s$quants[[2]][1,1],
      gamma_lo   = s$quants[[2]][2,1],
      gamma_hi   = s$quants[[2]][3,1],
      eta_est    = s$quants[[2]][1,2],
      q_est      = s$quants[[2]][1,3],
      z_elpd     = s$p,
      p_elpd     = pnorm(s$p)
    ))
  }

  X <- tryCatch(
    gwas_merge(exp_full, out_full,
               snp_name_cols = c("SNP","SNP"),
               beta_hat_cols = c("beta","beta"),
               se_cols       = c("se","se"),
               A1_cols       = c("effect_allele","effect_allele"),
               A2_cols       = c("other_allele","other_allele")),
    error = function(e){ cat("    merge fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(X) || nrow(X) < 1e5) {
    cat("    skip CAUSE (insufficient overlap, n=",
        ifelse(is.null(X),0,nrow(X)),")\n"); return(NULL)
  }

  # Clean extreme SNPs BEFORE nuisance estimation
  n_before <- nrow(X)
  X <- clean_cause_X(X)
  cat(sprintf("    SNPs after extreme-value cleaning: %d (was %d)\n",
              nrow(X), n_before))
  if (nrow(X) < 1e5) {
    cat("    too few SNPs after cleaning -> skip CAUSE\n"); return(NULL)
  }

  # ---- Estimate nuisance parameters (wall-clock timeout only) ----
  set.seed(123)
  vars <- sample(X$snp, size = min(2e5, nrow(X)), replace = FALSE)
  params <- tryCatch(
    withTimeout(
      est_cause_params(X, vars),
      cpu = Inf, elapsed = NUISANCE_TIMEOUT, onTimeout = "error"
    ),
    TimeoutException = function(e){
      cat("    [TIMEOUT] est_cause_params >", NUISANCE_TIMEOUT, "s -> skip\n"); NULL },
    error = function(e){
      cat("    nuisance fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(params)) return(NULL)

  # ---- ★ STRICT pval_thresh for PTSD / AUDIT_P / ANX ----
  is_strict <- (exp_name %in% STRICT_DISORDERS) || (out_name %in% STRICT_DISORDERS)
  pthresh   <- if (is_strict) 1e-4 else 1e-3
  cat(sprintf("    pval_thresh = %.0e  %s\n", pthresh,
              if (is_strict) "[strict-tightened]" else ""))

  # ---- Compute exposure p-values and pick top SNPs ----
  p1_vec  <- 2 * pnorm(-abs(X$beta_hat_1 / X$seb1))
  top_idx <- which(p1_vec < pthresh)
  cat("    candidate SNPs at p1 <", pthresh, ":", length(top_idx), "\n")

  if (length(top_idx) < 3) {
    relax <- pthresh * 10
    cat("    relax to p1 <", relax, "\n")
    top_idx <- which(p1_vec < relax)
    cat("    candidate SNPs at p1 <", relax, ":", length(top_idx), "\n")
  }
  if (length(top_idx) < 3) {
    cat("    still too few candidates -> skip CAUSE\n")
    return(NULL)
  }

  # ---- LD prune at r2<0.1 (per CAUSE authors) ----
  top_df <- tibble::tibble(rsid = X$snp[top_idx],
                           pval = p1_vec[top_idx])
  pruned <- tryCatch(
    ld_clump(top_df,
             clump_r2  = 0.1,
             clump_kb  = 10000,
             clump_p   = 1,
             plink_bin = PLINK, bfile = BFILE),
    error = function(e){ cat("    CAUSE prune fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(pruned) || nrow(pruned) < 3) {
    cat("    too few pruned SNPs (<3) -> skip CAUSE\n"); return(NULL)
  }
  cat("    SNPs after LD pruning (r2<0.1):", nrow(pruned), "\n")

  # ---- ★ Hard cap on posterior SNP count ----
  if (nrow(pruned) > POSTERIOR_SNP_CAP) {
    pruned <- pruned[order(pruned$pval), , drop = FALSE]
    pruned <- pruned[seq_len(POSTERIOR_SNP_CAP), , drop = FALSE]
    cat(sprintf("    capped posterior SNPs to top %d by p-value\n",
                POSTERIOR_SNP_CAP))
  }

  # ---- Run CAUSE (wall-clock timeout only) ----
  t0 <- Sys.time()
  res <- tryCatch(
    withTimeout(
      cause(X = X, variants = pruned$rsid, param_ests = params),
      cpu = Inf, elapsed = CAUSE_TIMEOUT, onTimeout = "error"
    ),
    TimeoutException = function(e){
      cat("    [TIMEOUT] cause() >", CAUSE_TIMEOUT, "s -> skip pair\n"); NULL },
    error = function(e){
      cat("    CAUSE fail:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) return(NULL)
  cat(sprintf("    cause() took %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  s <- summary(res, ci_size = 0.95)
  saveRDS(res, rds_path)

  tryCatch({
    pdf(file.path(outdir,"plots","cause", paste0(tag,"_cause.pdf")), w=8, h=6)
    print(plot(res, type="data")); print(plot(res)); dev.off()
  }, error = function(e) cat("    CAUSE plot fail:", conditionMessage(e), "\n"))

  data.frame(
    exposure   = exp_name, outcome = out_name,
    n_variants = nrow(pruned),
    gamma_est  = s$quants[[2]][1,1],
    gamma_lo   = s$quants[[2]][2,1],
    gamma_hi   = s$quants[[2]][3,1],
    eta_est    = s$quants[[2]][1,2],
    q_est      = s$quants[[2]][1,3],
    z_elpd     = s$p,
    p_elpd     = pnorm(s$p)
  )
}

# ============================================================================
# 6. 2SMR + MR-PRESSO
#    Resume: skip if results/{tag}_main.tsv already exists (non-empty)
#    If <3 SNPs after harmonise/steiger, retry with relaxed IVs (p<1e-5)
# ============================================================================
.do_2smr_inner <- function(exp_iv, out_full, exp_name, out_name){
  exp_fmt <- format_data(exp_iv, type="exposure",
            snp_col="SNP", beta_col="beta", se_col="se",
            effect_allele_col="effect_allele", other_allele_col="other_allele",
            eaf_col="eaf", pval_col="pval", samplesize_col="samplesize",
            phenotype_col="Phenotype")
  out_fmt <- format_data(out_full, type="outcome", snps=exp_fmt$SNP,
            snp_col="SNP", beta_col="beta", se_col="se",
            effect_allele_col="effect_allele", other_allele_col="other_allele",
            eaf_col="eaf", pval_col="pval", samplesize_col="samplesize",
            phenotype_col="Phenotype")
  if (nrow(out_fmt) < 3) {
    cat("   too few overlap (n_outcome=", nrow(out_fmt), ")\n"); return(NULL)
  }

  dat <- harmonise_data(exp_fmt, out_fmt, action=2)
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) < 3) { cat("   <3 SNPs after harmonise\n"); return(NULL) }

  dat <- tryCatch(steiger_filtering(dat),
                  error=function(e){ cat("   steiger fail:", conditionMessage(e),"\n"); dat })
  if ("steiger_dir" %in% names(dat))
    dat <- dat[dat$steiger_dir == TRUE | is.na(dat$steiger_dir), ]
  if (nrow(dat) < 3) { cat("   <3 SNPs after steiger\n"); return(NULL) }

  dat
}

run_2smr <- function(exp_iv, out_full, exp_name, out_name, outdir,
                     exp_full = NULL){
  tag <- paste0(exp_name, "_to_", out_name)
  cat(">>> 2SMR:", tag, "\n")

  # ---- RESUME ----
  main_path <- file.path(outdir, "results", paste0(tag,"_main.tsv"))
  if (file_done(main_path)) {
    cat("   [skip] already done ->", main_path, "\n")
    return(tryCatch(fread(main_path, data.table=FALSE),
                    error = function(e) NULL))
  }

  # ----- first attempt with given IVs (typically p<5e-8) -----
  dat <- NULL
  if (!is.null(exp_iv) && nrow(exp_iv) >= 3) {
    dat <- tryCatch(.do_2smr_inner(exp_iv, out_full, exp_name, out_name),
                    error = function(e){
                      cat("   first-pass error:", conditionMessage(e), "\n"); NULL })
  } else {
    cat("   no usable IV at default threshold\n")
  }

  # ----- retry with relaxed IVs (p<1e-5) if first attempt insufficient -----
  if (is.null(dat) && !is.null(exp_full)) {
    cat("   --> retrying with relaxed IV threshold p<1e-5\n")
    exp_iv2 <- get_iv(exp_full, p = 1e-5)
    if (!is.null(exp_iv2) && nrow(exp_iv2) >= 3 &&
        (is.null(exp_iv) || nrow(exp_iv2) > nrow(exp_iv))) {
      cat("   relaxed IV count:", nrow(exp_iv2), "\n")
      dat <- tryCatch(.do_2smr_inner(exp_iv2, out_full, exp_name, out_name),
                      error = function(e){
                        cat("   relaxed-pass error:", conditionMessage(e), "\n"); NULL })
      iv_relaxed_path <- file.path(outdir, "exposure_dat",
                                   paste0(exp_name,"_IVs_relaxed_1e5.tsv"))
      if (!file.exists(iv_relaxed_path))
        fwrite(exp_iv2, iv_relaxed_path, sep="\t")
    } else {
      cat("   relaxed IV extraction failed or did not increase set size\n")
    }
  }

  if (is.null(dat) || nrow(dat) < 3) {
    cat("   FINAL: <3 SNPs even after relaxation -> skip pair\n"); return(NULL)
  }

  fs <- f_stat(dat)
  cat(sprintf("   k=%d  R2=%s  F=%s\n",
              fs$k, fmt_num(fs$R2, 4), fmt_num(fs$F, 2)))

  res <- mr(dat, method_list=c("mr_ivw","mr_egger_regression",
                               "mr_weighted_median","mr_weighted_mode",
                               "mr_ivw_mre"))
  res$exposure <- exp_name; res$outcome <- out_name
  res$F_stat <- fs$F; res$R2 <- fs$R2; res$k_SNP <- fs$k

  het   <- mr_heterogeneity(dat)
  pleio <- mr_pleiotropy_test(dat)
  loo   <- mr_leaveoneout(dat)

  presso <- NULL
  if (nrow(dat) >= 4) {
    presso <- tryCatch(
      run_mr_presso(dat, NbDistribution=1000, SignifThreshold=0.05),
      error=function(e){ cat("   PRESSO fail:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(presso)) cat("   PRESSO done\n")
  } else {
    cat("   PRESSO skipped: need >=4 SNPs (have", nrow(dat), ")\n")
  }

  fwrite(res,   main_path,                                              sep="\t")
  fwrite(het,   file.path(outdir,"sensitivity", paste0(tag,"_het.tsv")),  sep="\t")
  fwrite(pleio, file.path(outdir,"sensitivity", paste0(tag,"_pleio.tsv")),sep="\t")
  fwrite(loo,   file.path(outdir,"sensitivity", paste0(tag,"_loo.tsv")),  sep="\t")

  if (!is.null(presso) && length(presso) > 0) {
    pres <- presso[[1]]

    pr_main <- pres$`Main MR results`
    if (!is.null(pr_main)) {
      pr_main$exposure <- exp_name
      pr_main$outcome  <- out_name
      fwrite(pr_main, file.path(outdir,"sensitivity",
                                paste0(tag,"_presso_main.tsv")), sep="\t")
    }

    gt <- pres$`MR-PRESSO results`$`Global Test`
    if (!is.null(gt)) {
      gt_df <- data.frame(exposure=exp_name, outcome=out_name,
                          RSSobs = ifelse(is.null(gt$RSSobs), NA, gt$RSSobs),
                          Pvalue = ifelse(is.null(gt$Pvalue), NA, gt$Pvalue))
      fwrite(gt_df, file.path(outdir,"sensitivity",
                              paste0(tag,"_presso_global.tsv")), sep="\t")
    }

    ot <- pres$`MR-PRESSO results`$`Outlier Test`
    if (!is.null(ot) && is.data.frame(ot) && nrow(ot) > 0) {
      ot$SNP      <- dat$SNP
      ot$exposure <- exp_name
      ot$outcome  <- out_name
      fwrite(ot, file.path(outdir,"sensitivity",
                           paste0(tag,"_presso_outlier.tsv")), sep="\t")
    }

    dist_t <- pres$`MR-PRESSO results`$`Distortion Test`
    if (!is.null(dist_t) && length(dist_t) > 0) {
      oi <- dist_t$`Outliers Indices`
      n_out <- if (is.null(oi) || identical(as.character(oi)[1],"No significant outliers")) 0 else length(oi)
      dist_df <- data.frame(
        exposure         = exp_name, outcome = out_name,
        n_outliers       = n_out,
        outliers_indices = ifelse(n_out==0, NA, paste(as.character(oi), collapse=",")),
        distortion_coef  = ifelse(is.null(dist_t$`Distortion Coefficient`), NA, dist_t$`Distortion Coefficient`),
        Pvalue           = ifelse(is.null(dist_t$Pvalue), NA, dist_t$Pvalue)
      )
      fwrite(dist_df, file.path(outdir,"sensitivity",
                                paste0(tag,"_presso_distortion.tsv")), sep="\t")
    }

    saveRDS(presso, file.path(outdir,"sensitivity", paste0(tag,"_presso.rds")))
  }

  saveRDS(dat, file.path(outdir,"harmonised", paste0(tag,".rds")))

  tryCatch({
    p1 <- mr_scatter_plot(res, dat)[[1]]
    p2 <- mr_funnel_plot(mr_singlesnp(dat))[[1]]
    p3 <- mr_forest_plot(mr_singlesnp(dat))[[1]]
    ggsave(file.path(outdir,"plots",paste0(tag,"_scatter.pdf")), p1, w=6, h=5)
    ggsave(file.path(outdir,"plots",paste0(tag,"_funnel.pdf")),  p2, w=6, h=5)
    ggsave(file.path(outdir,"plots",paste0(tag,"_forest.pdf")),  p3, w=6, h=8)
  }, error=function(e) cat("   plot fail:", conditionMessage(e), "\n"))

  res
}

# Helper: check if both 2SMR and CAUSE for a pair are fully done
pair_fully_done <- function(exp_name, out_name, outdir){
  tag <- paste0(exp_name, "_to_", out_name)
  main_done  <- file_done(file.path(outdir,"results", paste0(tag,"_main.tsv")))
  cause_done <- file_done(file.path(outdir,"results", paste0(tag,"_cause.rds")))
  main_done && cause_done
}

# Helper: load cached CAUSE summary row (used when CAUSE already done)
load_cause_summary <- function(exp_name, out_name, outdir){
  tag <- paste0(exp_name, "_to_", out_name)
  rds_path <- file.path(outdir, "results", paste0(tag,"_cause.rds"))
  if (!file_done(rds_path)) return(NULL)
  res <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  s <- tryCatch(summary(res, ci_size = 0.95), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  data.frame(
    exposure   = exp_name, outcome = out_name,
    n_variants = NA_integer_,
    gamma_est  = s$quants[[2]][1,1],
    gamma_lo   = s$quants[[2]][2,1],
    gamma_hi   = s$quants[[2]][3,1],
    eta_est    = s$quants[[2]][1,2],
    q_est      = s$quants[[2]][1,3],
    z_elpd     = s$p,
    p_elpd     = pnorm(s$p)
  )
}

# ============================================================================
# 7. Main loop
# ============================================================================
cat("\n########## STEP 1. Read all GWAS ##########\n")
gwas <- list()
for (nm in names(mvpa_list)) gwas[[nm]] <- read_gwas(mvpa_list[[nm]], nm, file.path(opt$out,"cache"))
for (nm in names(disorders)) gwas[[nm]] <- read_gwas(disorders[[nm]], nm, file.path(opt$out,"cache"))

cat("\n########## STEP 2. Extract instruments ##########\n")
iv <- list()
for (nm in names(gwas)){
  cat("  ", nm, "\n")
  iv_path <- file.path(opt$out,"exposure_dat", paste0(nm,"_IVs.tsv"))
  iv[[nm]] <- get_iv_with_fallback(gwas[[nm]], nm, iv_path)
  if (is.null(iv[[nm]])) cat("    no IVs at any threshold\n")
}

# ----- FORWARD : Psychiatric disorders -> MVPA -----
cat("\n########## STEP 3. FORWARD  MR + CAUSE : Disorders -> MVPA ##########\n")
fwd_mr <- list(); fwd_cause <- list()
for (dz in names(disorders)){
  for (mv in names(mvpa_list)){
    if (is.null(gwas[[mv]])) next

    # Fast-path: pair fully done -> just load summaries, no heavy IO
    if (pair_fully_done(dz, mv, opt$out)) {
      cat(sprintf("[skip-pair] %s -> %s already complete\n", dz, mv))
      mp <- file.path(opt$out,"results", paste0(dz,"_to_",mv,"_main.tsv"))
      r  <- tryCatch(fread(mp, data.table=FALSE), error = function(e) NULL)
      if (!is.null(r)) fwd_mr[[paste0(dz,"_",mv)]] <- r
      cc <- load_cause_summary(dz, mv, opt$out)
      if (!is.null(cc)) fwd_cause[[paste0(dz,"_",mv)]] <- cc
      next
    }

    # 2SMR
    r <- tryCatch(
      run_2smr(iv[[dz]], gwas[[mv]], dz, mv, opt$out, exp_full = gwas[[dz]]),
      error=function(e){cat("   ERR:",conditionMessage(e),"\n");NULL})
    if (!is.null(r)) fwd_mr[[paste0(dz,"_",mv)]] <- r

    # CAUSE
    cc <- tryCatch(run_cause(gwas[[dz]], gwas[[mv]], dz, mv, opt$out),
                   error=function(e){cat("   CAUSE ERR:",conditionMessage(e),"\n");NULL})
    if (!is.null(cc)) fwd_cause[[paste0(dz,"_",mv)]] <- cc
  }
}

# ----- REVERSE : MVPA -> Psychiatric disorders -----
cat("\n########## STEP 4. REVERSE  MR + CAUSE : MVPA -> Disorders ##########\n")
rev_mr <- list(); rev_cause <- list()
for (mv in names(mvpa_list)){
  for (dz in names(disorders)){
    if (is.null(gwas[[dz]])) next

    if (pair_fully_done(mv, dz, opt$out)) {
      cat(sprintf("[skip-pair] %s -> %s already complete\n", mv, dz))
      mp <- file.path(opt$out,"results", paste0(mv,"_to_",dz,"_main.tsv"))
      r  <- tryCatch(fread(mp, data.table=FALSE), error = function(e) NULL)
      if (!is.null(r)) rev_mr[[paste0(mv,"_",dz)]] <- r
      cc <- load_cause_summary(mv, dz, opt$out)
      if (!is.null(cc)) rev_cause[[paste0(mv,"_",dz)]] <- cc
      next
    }

    r <- tryCatch(
      run_2smr(iv[[mv]], gwas[[dz]], mv, dz, opt$out, exp_full = gwas[[mv]]),
      error=function(e){cat("   ERR:",conditionMessage(e),"\n");NULL})
    if (!is.null(r)) rev_mr[[paste0(mv,"_",dz)]] <- r

    cc <- tryCatch(run_cause(gwas[[mv]], gwas[[dz]], mv, dz, opt$out),
                   error=function(e){cat("   CAUSE ERR:",conditionMessage(e),"\n");NULL})
    if (!is.null(cc)) rev_cause[[paste0(mv,"_",dz)]] <- cc
  }
}

# ============================================================================
# 8. Summary
# ============================================================================
cat("\n########## STEP 5. Summarise ##########\n")
summarise <- function(lst, direction){
  if (length(lst)==0) return(NULL)
  out <- do.call(rbind, lst)
  out$direction <- direction
  out
}

mr_parts <- list()
fwd_sum <- summarise(fwd_mr, "forward")
rev_sum <- summarise(rev_mr, "reverse")
if (!is.null(fwd_sum)) mr_parts[["forward"]] <- fwd_sum
if (!is.null(rev_sum)) mr_parts[["reverse"]] <- rev_sum
if (length(mr_parts) > 0) {
  all_mr <- do.call(rbind, mr_parts)
  fwrite(all_mr, file.path(opt$out,"results","MR_summary_all.tsv"), sep="\t")
}

cause_parts <- list()
if (length(fwd_cause) > 0)
  cause_parts[["forward"]] <- cbind(do.call(rbind, fwd_cause), direction="forward")
if (length(rev_cause) > 0)
  cause_parts[["reverse"]] <- cbind(do.call(rbind, rev_cause), direction="reverse")
if (length(cause_parts) > 0) {
  cs <- do.call(rbind, cause_parts)
  fwrite(cs, file.path(opt$out,"results","CAUSE_summary_all.tsv"), sep="\t")
}

cat("\nAll done. See: ", file.path(opt$out,"results"), "\n")

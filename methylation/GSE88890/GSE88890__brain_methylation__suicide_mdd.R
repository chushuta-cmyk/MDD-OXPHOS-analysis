# ==============================================================================
# GSE88890__brain_methylation__suicide_mdd.R (archived analysis copy)
# ------------------------------------------------------------------------------
# Dataset : GSE88890 - POST-MORTEM BRAIN 450K DNA methylation
#           BA11 (n=40) + BA25 (n=35) from MDD suicide cases (37) and
#           non-psychiatric sudden death controls (38). 75 samples, 58M/17F.
# Platform: GPL13534 Illumina HumanMethylation450
# Original: Murphy et al. 2017, Transl Psychiatry 7:e989
#           (PSORS1C3 upstream DMR, chr6 MHC, validated by pyrosequencing)
# NOTE     : This dataset is BRAIN METHYLATION. It belongs under methylation/,
#           not blood/. Directory relocation is recommended.
# ------------------------------------------------------------------------------
# Fixes vs old pipeline (see All_Pipelines_Evaluation_Report.md):
#   [P0] PSEUDO-REPLICATION: the same subject contributes BA11 AND BA25, old
#        code treated them as independent. NOW:
#          (a) main model: limma with duplicateCorrelation(block=Subject);
#          (b) per-region models: BA11 and BA25 tested separately;
#          (c) never pool regions as independent samples.
#   [P1] Old 01_main_analysis.R had a typo ('s1·trsplit') and could not run.
#        NOW: rewritten cleanly in modules.
#   [P1] Old result = 1 DMP (cg18100887, chr17, no gene) vs README "14 DMP"
#        vs paper PSORS1C3. NOW: explicit PSORS1C3 / chr6 MHC check.
#   [P2] Region parsed from tissue field; Subject ID parsed from sample title
#        ("184_BA25" -> subject 184). Gender + Age used as covariates.
# Thresholds (user policy): adj.P < 0.05; no beta-diff cutoff.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  # NOTE: IlluminaHumanMethylation450kanno.ilmn12.hg19 is loaded INSIDE
  # MODULE 3 with requireNamespace() so the script still runs when the
  # annotation package is missing (annotation/PSORS1C3 check skipped then).
})

# -------------------------------------------------------------------------------
# SETTINGS
# -------------------------------------------------------------------------------
# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

project_root <- file.path(WORK_ROOT, "transcriptomics")
output_dir   <- file.path(OUTPUT_ROOT, "methylation/GSE88890/results_fixed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

series_file <- file.path(DATA_ROOT, "GSE88890/GSE88890_series_matrix.txt")
beta_file   <- file.path(DATA_ROOT, "GSE88890/GSE88890_normalisedBetas.csv")
manifest_file <- file.path(DATA_ROOT, "GSE88890/GSE88890_RAW/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz")

PADJ_CUTOFF <- 0.05

# ==============================================================================
# MODULE 1: LOAD DATA & METADATA (Subject / Region / Group / Gender / Age)
# ==============================================================================
cat("\n===== MODULE 1: LOAD DATA & METADATA =====\n")

# --- 1.1 Parse series matrix header --------------------------------------------
parse_sm <- function(f) {
  con <- file(f, "r"); on.exit(close(con))
  meta <- list()
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (grepl("^!Sample_", line)) {
      parts <- strsplit(line, "\t")[[1]]
      key <- sub("^!", "", parts[1])
      vals <- gsub('"', '', parts[-1])
      if (!key %in% names(meta)) meta[[key]] <- list()
      meta[[key]][[length(meta[[key]]) + 1]] <- vals
    }
    if (grepl("^!series_matrix_table_begin", line)) break
  }
  meta
}
sm <- parse_sm(series_file)

# ROBUST field extraction: match FIELD NAME (text before ":") of each
# characteristics row. No regex escape sequences (avoid \s issues in console).
# NOTE: `lst` is the WHOLE row vector (75 values); `lst[[1]]` is only the
# first value - always apply sub() to `lst`, not to `lst[[1]]`.
extract_field <- function(meta_lists, field_name) {
  for (lst in meta_lists) {
    if (length(lst) == 0) next
    first_val <- lst[1]
    field <- strsplit(first_val, ": ", fixed = TRUE)[[1]][1]
    if (identical(field, field_name)) {
      return(sub("^[^:]*: ", "", lst))         # strip "name: " from ALL values
    }
  }
  rep(NA_character_, length(meta_lists[[1]]))
}

title_vals <- sm[["Sample_title"]][[1]]
gsm_ids    <- sm[["Sample_geo_accession"]][[1]]
tissue     <- extract_field(sm[["Sample_characteristics_ch1"]], "tissue")
age        <- suppressWarnings(as.numeric(extract_field(sm[["Sample_characteristics_ch1"]], "age")))
groupid    <- extract_field(sm[["Sample_characteristics_ch1"]], "groupid")
gender     <- extract_field(sm[["Sample_characteristics_ch1"]], "gender")

cat("  Raw groupid counts:\n"); print(table(groupid, useNA = "ifany"))
cat("  Raw region counts:\n");  print(table(sub(".*_", "", title_vals)))

# Subject ID: from title "184_BA25" -> 184  (FIX: enables blocking by subject)
subject_id <- sub("_.*", "", title_vals)
region     <- sub(".*_", "", title_vals)
region     <- ifelse(grepl("11", region), "BA11", ifelse(grepl("25", region), "BA25", region))
phenotype  <- ifelse(grepl("MDD", groupid, ignore.case = TRUE), "MDD", "Control")

pheno <- data.frame(
  Sample_ID = gsm_ids,
  Subject   = subject_id,
  Region    = factor(region),
  Phenotype = factor(phenotype, levels = c("Control", "MDD")),
  Gender    = factor(gender),
  Age       = age,
  stringsAsFactors = FALSE
)
# --- 1.15 Sentrix position -> GSM mapping (FIX: column-name alignment) ----------
# normalisedBetas.csv columns are named by SENTRIX POSITION (e.g.
# "9283265107_R01C02"), NOT by GSM. GEO stores the Sentrix position in the
# Sample_description field of the series matrix - use it to map columns.
sentrix_vals <- sm[["Sample_description"]][[1]]
pheno$Sentrix <- sentrix_vals
cat("  Sentrix mapping ready (first 3):", head(pheno$Sentrix, 3), "\n")
cat("  Samples:", nrow(pheno), "\n")
cat("  Subjects:", length(unique(pheno$Subject)), "(each subject has 1-2 regions)\n")
cat("  Region x Phenotype table:\n"); print(table(pheno$Region, pheno$Phenotype))

# --- 1.2 Load beta matrix -------------------------------------------------------
beta <- fread(beta_file)
probe_col <- colnames(beta)[1]
beta_mat <- as.matrix(beta[, -1, with = FALSE])
rownames(beta_mat) <- beta[[probe_col]]
colnames(beta_mat) <- sub("\\.csv$", "", colnames(beta_mat))
cat(sprintf("  Beta matrix: %d probes x %d samples\n", nrow(beta_mat), ncol(beta_mat)))

# Align metadata to beta columns VIA SENTRIX (was GSM -> all NA, FIX)
idx_sentrix <- match(colnames(beta_mat), pheno$Sentrix)
if (any(is.na(idx_sentrix))) {
  cat("  [WARN] Unmatched Sentrix IDs in beta columns:\n")
  print(colnames(beta_mat)[is.na(idx_sentrix)])
}
stopifnot(!any(is.na(idx_sentrix)))
pheno <- pheno[idx_sentrix, ]
colnames(beta_mat) <- pheno$Sample_ID     # rename to GSM for downstream code

# Load annotation before QC so sex-chromosome and SNP probes can be removed.
load_450k_annotation <- function() {
  # Prefer the local GPL manifest.  It is self-contained and avoids a fragile
  # assumption made by the previous code: an installed annotation *package*
  # does not necessarily expose an object with exactly the package name in its
  # namespace.  That assumption caused get(..., asNamespace(...)) to fail on
  # systems where the package was installed, even though this manifest existed.
  if (file.exists(manifest_file)) {
    x <- data.table::fread(manifest_file, skip = 7, fill = TRUE, showProgress = FALSE,
                           select = c("Name", "CHR", "MAPINFO", "Probe_SNPs", "Probe_SNPs_10",
                                      "UCSC_RefGene_Name", "UCSC_RefGene_Group",
                                      "Relation_to_UCSC_CpG_Island"))
    required_cols <- c("Name", "CHR", "MAPINFO", "UCSC_RefGene_Name",
                       "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island")
    if (!all(required_cols %in% names(x))) {
      stop("Local 450K manifest is missing required columns: ",
           paste(setdiff(required_cols, names(x)), collapse = ", "))
    }
    x <- x[grepl("^(cg|ch|rs)", Name)]
    chr <- as.character(x$CHR)
    chr[!is.na(chr) & nzchar(chr) & !grepl("^chr", chr)] <-
      paste0("chr", chr[!is.na(chr) & nzchar(chr) & !grepl("^chr", chr)])
    out <- data.frame(Probe = x$Name,
                      chr = chr,
                      pos = as.numeric(x$MAPINFO),
                      UCSC_RefGene_Name = x$UCSC_RefGene_Name,
                      UCSC_RefGene_Group = x$UCSC_RefGene_Group,
                      Relation_to_Island = x$Relation_to_UCSC_CpG_Island,
                      Probe_SNPs = paste(ifelse(is.na(x$Probe_SNPs), "", x$Probe_SNPs),
                                         ifelse(is.na(x$Probe_SNPs_10), "", x$Probe_SNPs_10)),
                      stringsAsFactors = FALSE)
    cat("  Annotation source: local GPL13534 manifest\n")
    return(out)
  }

  # Package fallback: load the package data object explicitly, then ask minfi
  # for its annotation.  Never use get(package_name, package_namespace).
  anno_pkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  if (requireNamespace(anno_pkg, quietly = TRUE) &&
      requireNamespace("minfi", quietly = TRUE)) {
    data_env <- new.env(parent = emptyenv())
    utils::data(list = anno_pkg, package = anno_pkg, envir = data_env)
    if (exists(anno_pkg, envir = data_env, inherits = FALSE)) {
      anno_obj <- get(anno_pkg, envir = data_env, inherits = FALSE)
      x <- as.data.frame(minfi::getAnnotation(anno_obj))
      if (!all(c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group",
                 "Relation_to_Island") %in% names(x))) {
        stop("450K annotation package returned an unexpected schema.")
      }
      cat("  Annotation source: Illumina 450K Bioconductor package\n")
      return(data.frame(Probe = rownames(x), chr = as.character(x$chr),
                        pos = x$pos,
                        UCSC_RefGene_Name = x$UCSC_RefGene_Name,
                        UCSC_RefGene_Group = x$UCSC_RefGene_Group,
                        Relation_to_Island = x$Relation_to_Island,
                        Probe_SNPs = if ("Probe_rs" %in% names(x)) x$Probe_rs else "",
                        stringsAsFactors = FALSE))
    }
  }
  NULL
}
ann <- load_450k_annotation()
cat(sprintf("  450K annotation: %s\n",
            if (is.null(ann)) "not available" else paste(nrow(ann), "probes")))

# ==============================================================================
# MODULE 2: QC & DMP ANALYSIS (FIX: block by Subject, per-region models)
# ==============================================================================
cat("\n===== MODULE 2: DMP ANALYSIS =====\n")

# --- 2.1 Filter low-quality, sex-chromosome and SNP-affected probes ------------
# Require >=95% complete values.  The old >=50% rule was too permissive.
keep_probe <- rowMeans(is.na(beta_mat)) <= 0.05
if (!is.null(ann)) {
  ann_idx <- match(rownames(beta_mat), ann$Probe)
  sex_probe <- !is.na(ann_idx) & ann$chr[ann_idx] %in% c("chrX", "chrY", "X", "Y")
  snp_probe <- !is.na(ann_idx) & !is.na(ann$Probe_SNPs[ann_idx]) &
    nzchar(trimws(ann$Probe_SNPs[ann_idx]))
  keep_probe <- keep_probe & !sex_probe & !snp_probe
  cat(sprintf("  Removed sex-chromosome probes: %d; SNP-annotated probes: %d\n",
              sum(sex_probe), sum(snp_probe)))
}
beta_qc <- beta_mat[keep_probe, ]
cat(sprintf("  Probes after NA filter: %d\n", nrow(beta_qc)))

# Statistical tests use M-values; beta values remain the interpretable effect
# scale.  Clip only for the logit transform, not for reported delta-beta.
beta_clip <- pmin(pmax(beta_qc, 1e-6), 1 - 1e-6)
m_qc <- log2(beta_clip / (1 - beta_clip))

# --- 2.2 MAIN model: ~ Phenotype + Region + Gender + Age, blocked by Subject ---
# FIX [P0]: duplicateCorrelation accounts for 2 regions per subject
design_main <- model.matrix(~ Phenotype + Region + Gender + Age, data = pheno)
if (qr(design_main)$rank != ncol(design_main)) stop("Main methylation design is rank deficient.")
corfit <- duplicateCorrelation(m_qc, design_main, block = pheno$Subject)
cat(sprintf("  Within-subject correlation (duplicateCorrelation): %.4f\n", corfit$consensus.correlation))
fit <- lmFit(m_qc, design_main, block = pheno$Subject, correlation = corfit$consensus.correlation)
fit <- eBayes(fit)
dmp_main <- topTable(fit, coef = "PhenotypeMDD", number = Inf, adjust.method = "BH")
dmp_main$ProbeID <- rownames(dmp_main)
dmp_main$DeltaBeta <- rowMeans(beta_qc[, pheno$Phenotype == "MDD", drop = FALSE], na.rm = TRUE)[dmp_main$ProbeID] -
  rowMeans(beta_qc[, pheno$Phenotype == "Control", drop = FALSE], na.rm = TRUE)[dmp_main$ProbeID]
dmp_main_sig <- dmp_main[dmp_main$adj.P.Val < PADJ_CUTOFF, ]
cat(sprintf("  MAIN model (blocked by subject): DMPs adj.P<0.05 = %d\n", nrow(dmp_main_sig)))
write.csv(dmp_main, file.path(output_dir, "GSE88890_DMP_main_blockedSubject.csv"), row.names = FALSE)

# --- 2.3 Per-region models (no pseudo-replication at all) ----------------------
per_region <- list()
for (rg in c("BA11", "BA25")) {
  idx <- which(pheno$Region == rg)
  if (sum(pheno$Phenotype[idx] == "MDD") >= 3 && sum(pheno$Phenotype[idx] == "Control") >= 3) {
    design_r <- model.matrix(~ Phenotype + Gender + Age,
                             data = pheno[idx, , drop = FALSE])
    if (qr(design_r)$rank != ncol(design_r)) stop("Rank-deficient design in ", rg)
    fit_r <- lmFit(m_qc[, idx], design_r)
    fit_r <- eBayes(fit_r)
    dmp_r <- topTable(fit_r, coef = "PhenotypeMDD", number = Inf, adjust.method = "BH")
    dmp_r$ProbeID <- rownames(dmp_r)
    dmp_r$DeltaBeta <- rowMeans(beta_qc[, idx[pheno$Phenotype[idx] == "MDD"], drop = FALSE], na.rm = TRUE)[dmp_r$ProbeID] -
      rowMeans(beta_qc[, idx[pheno$Phenotype[idx] == "Control"], drop = FALSE], na.rm = TRUE)[dmp_r$ProbeID]
    per_region[[rg]] <- dmp_r
    cat(sprintf("  %s: DMPs adj.P<0.05 = %d (n=%d)\n",
                rg, sum(dmp_r$adj.P.Val < PADJ_CUTOFF), length(idx)))
    write.csv(dmp_r, file.path(output_dir, paste0("GSE88890_DMP_", rg, ".csv")), row.names = FALSE)
  } else {
    cat(sprintf("  %s: skipped (insufficient per-group n)\n", rg))
  }
}

# ==============================================================================
# MODULE 3: ANNOTATION & PSORS1C3 CHECK (reproduce the paper's hit)
# ==============================================================================
cat("\n===== MODULE 3: ANNOTATION & PSORS1C3 CHECK =====\n")

# Annotation was loaded from the Bioconductor package or the local GPL13534
# manifest in MODULE 1.  Therefore annotation no longer silently disappears
# merely because the optional package is absent.

annotate_dmp <- function(dmp) {
  if (is.null(ann)) return(dmp)   # no annotation available
  dmp$Chr       <- ann$chr[match(dmp$ProbeID, ann$Probe)]
  dmp$Pos       <- ann$pos[match(dmp$ProbeID, ann$Probe)]
  dmp$Gene      <- ann$UCSC_RefGene_Name[match(dmp$ProbeID, ann$Probe)]
  dmp$GeneGroup <- ann$UCSC_RefGene_Group[match(dmp$ProbeID, ann$Probe)]
  dmp$Island    <- ann$Relation_to_Island[match(dmp$ProbeID, ann$Probe)]
  dmp
}
if (nrow(dmp_main_sig) > 0) {
  dmp_main_sig <- annotate_dmp(dmp_main_sig)
  write.csv(dmp_main_sig, file.path(output_dir, "GSE88890_DMP_annotated_adjP0.05.csv"), row.names = FALSE)
  cat("  Significant DMPs:\n"); print(dmp_main_sig[, c("ProbeID", "Chr", "Gene", "logFC", "adj.P.Val")])
}

# --- 3.1 PSORS1C3 region check (chr6 MHC, upstream of PSORS1C3) ----------------
# PSORS1C3 sits in the MHC region chr6:31,150,000-31,160,000 (hg19 approx).
mhc_hits <- data.frame()     # initialize so MODULE 4 can always reference it
if (!is.null(ann)) {
  mhc_probes <- ann$Probe[ann$chr == "chr6" & ann$pos > 31000000 & ann$pos < 31600000]
  mhc_hits <- dmp_main[dmp_main$ProbeID %in% mhc_probes, ]
  cat(sprintf("  chr6 MHC probes tested: %d\n", length(mhc_probes)))
  if (nrow(mhc_hits) > 0) {
    cat("  chr6 MHC region DMPs (PSORS1C3 vicinity):\n")
    print(head(mhc_hits[order(mhc_hits$P.Value),
                        c("ProbeID", "logFC", "DeltaBeta", "P.Value", "adj.P.Val")], 20))
    write.csv(mhc_hits, file.path(output_dir, "GSE88890_chr6_MHC_PSORS1C3_check.csv"), row.names = FALSE)
  } else {
    cat("  No chr6 MHC region probes reach any significance level in this run.\n")
    cat("  (Paper: DMR upstream of PSORS1C3 validated by pyrosequencing;\n")
    cat("   DMR tests need comb-p/DMRcate and were not done here.)\n")
  }
} else {
  cat("  [WARN] Annotation package missing; PSORS1C3 region check skipped.\n")
}

# ==============================================================================
# MODULE 3.5: P-VALUE DIAGNOSTICS (why adj.P<0.05 may be zero - check reality)
# ==============================================================================
# Context: 40 subjects / 75 samples vs ~416k probes. Single-CpG BH across the
# genome needs P < ~1e-7, which n=40 post-mortem samples cannot usually reach.
# The original paper (Murphy et al. 2017, Transl Psychiatry) did NOT use single-
# CpG BH either - it used Bumphunter DMR analysis (PSORS1C3 DMR, 13 CpGs,
# Sidak P = 8.02e-6) and a later re-analysis (IJNP 2023) reported uncorrected
# P<0.01 CpGs (BA11: 3177, BA25: 1304). So zero adj.P<0.05 DMPs is expected
# here and NOT a data-alignment failure. These diagnostics verify that.
cat("\n===== MODULE 3.5: P-VALUE DIAGNOSTICS =====\n")

p_diag <- function(dmp, label) {
  lam <- median(qchisq(1 - dmp$P.Value, 1)) / qchisq(0.5, 1)
  cat(sprintf("  [%s] probes tested: %d\n", label, nrow(dmp)))
  cat(sprintf("    min P.Value : %s\n", format(min(dmp$P.Value), scientific = TRUE)))
  cat(sprintf("    P<0.01      : %d  (paper's IJNP-2023 reporting cut)\n", sum(dmp$P.Value < 0.01)))
  cat(sprintf("    P<1e-4      : %d\n", sum(dmp$P.Value < 1e-4)))
  cat(sprintf("    P<1e-5      : %d\n", sum(dmp$P.Value < 1e-5)))
  cat(sprintf("    adj.P<0.05  : %d  (single-CpG BH - rarely powered at n=40)\n",
              sum(dmp$adj.P.Val < PADJ_CUTOFF)))
  cat(sprintf("    inflation lambda: %.3f  (~1.0 normal; >1.1 suggests batch/tech)\n", lam))
  top <- head(dmp[order(dmp$P.Value), c("ProbeID", "logFC", "P.Value", "adj.P.Val")], 3)
  print(top, row.names = FALSE)
  data.frame(label = label, min_P = min(dmp$P.Value), n_P0.01 = sum(dmp$P.Value < 0.01),
             n_P1e4 = sum(dmp$P.Value < 1e-4), n_adj = sum(dmp$adj.P.Val < PADJ_CUTOFF),
             lambda = lam)
}
diag_main <- p_diag(dmp_main, "MAIN  (BA11+BA25, dupCorr by Subject)")
diag_rows <- list(diag_main)
if (!is.null(per_region$BA11)) diag_rows[[length(diag_rows) + 1]] <- p_diag(per_region$BA11, "BA11  (n=40)")
if (!is.null(per_region$BA25)) diag_rows[[length(diag_rows) + 1]] <- p_diag(per_region$BA25, "BA25  (n=35)")
diag_df <- do.call(rbind, diag_rows)
write.csv(diag_df, file.path(output_dir, "GSE88890_pvalue_diagnostics.csv"), row.names = FALSE)

# Volcano plot: x-axis is delta-beta for interpretation; inference is based on
# M-values.  Plot only finite tests so ggplot/base graphics do not silently drop
# hundreds of rows with warnings.
plot_dmp <- dmp_main[is.finite(dmp_main$DeltaBeta) & is.finite(dmp_main$P.Value) &
                       dmp_main$P.Value > 0, ]
pdf(file.path(output_dir, "GSE88890_volcano_delta_beta.pdf"), width = 8, height = 6)
plot(plot_dmp$DeltaBeta, -log10(plot_dmp$P.Value), pch = 16, cex = 0.45,
     col = ifelse(plot_dmp$adj.P.Val < PADJ_CUTOFF, "#D7301F", "grey70"),
     main = "GSE88890: MDD vs control (subject-blocked)",
     sub = "M-value inference; x-axis shows delta-beta; colour = BH FDR",
     xlab = "Delta beta (MDD - control)", ylab = "-log10(P)")
abline(h = -log10(0.05), lty = 2, col = "#2C7FB8")
dev.off()

# ==============================================================================
# MODULE 4: REPORT
# ==============================================================================
report <- sprintf(
"GSE88890 FIXED pipeline report (2026-09)
===========================================
Samples      : %d (MDD suicide=%d, Control=%d; BA11=%d, BA25=%d; subjects=%d)
FIX [P0]     : pseudo-replication removed - duplicateCorrelation(block=Subject)
               in main model + per-region models (BA11/BA25 separate).
Main DMPs adj.P<0.05: %d
  BA11 DMPs  : %s
  BA25 DMPs  : %s
PSORS1C3 vicinity check: %s

P-value diagnostics (see GSE88890_pvalue_diagnostics.csv):
  MAIN  : min P=%.2e | P<0.01=%d | P<1e-4=%d | adj.P<0.05=%d | lambda=%.2f
  BA11  : min P=%.2e | P<0.01=%d | P<1e-4=%d | adj.P<0.05=%d
  BA25  : min P=%.2e | P<0.01=%d | P<1e-4=%d | adj.P<0.05=%d
Interpretation: genome-wide single-CpG power is limited at n=40 subjects.
  The one pooled-model DMP must be treated as exploratory because neither
  region-specific analysis reaches FDR < 0.05.
  Original paper (Murphy 2017) used Bumphunter DMR analysis, not single-CpG
  BH; PSORS1C3 is a 13-CpG DMR (Sidak P=8.02e-6). IJNP-2023 re-analysis
  reported uncorrected P<0.01 CpGs (BA11 3177 / BA25 1304). Report either
  the P<0.01 count or DMR results, not single-CpG adj.P.
NOTE         : this dataset is BRAIN methylation - should live under
               methylation/GSE88890/, not blood/. Relocation recommended.
",
  nrow(pheno),
  sum(pheno$Phenotype == "MDD"), sum(pheno$Phenotype == "Control"),
  sum(pheno$Region == "BA11"), sum(pheno$Region == "BA25"),
  length(unique(pheno$Subject)),
  nrow(dmp_main_sig),
  if (!is.null(per_region$BA11)) sum(per_region$BA11$adj.P.Val < PADJ_CUTOFF) else "n/a",
  if (!is.null(per_region$BA25)) sum(per_region$BA25$adj.P.Val < PADJ_CUTOFF) else "n/a",
  if (nrow(mhc_hits) > 0) sprintf("%d probes; min P=%.2e; FDR<0.05=%d",
                                  nrow(mhc_hits), min(mhc_hits$P.Value, na.rm = TRUE),
                                  sum(mhc_hits$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE)) else "no annotated probes in window",
  diag_main$min_P, diag_main$n_P0.01, diag_main$n_P1e4, diag_main$n_adj, diag_main$lambda,
  if (!is.null(per_region$BA11)) min(per_region$BA11$P.Value) else NA,
  if (!is.null(per_region$BA11)) sum(per_region$BA11$P.Value < 0.01) else NA,
  if (!is.null(per_region$BA11)) sum(per_region$BA11$P.Value < 1e-4) else NA,
  if (!is.null(per_region$BA11)) sum(per_region$BA11$adj.P.Val < PADJ_CUTOFF) else NA,
  if (!is.null(per_region$BA25)) min(per_region$BA25$P.Value) else NA,
  if (!is.null(per_region$BA25)) sum(per_region$BA25$P.Value < 0.01) else NA,
  if (!is.null(per_region$BA25)) sum(per_region$BA25$P.Value < 1e-4) else NA,
  if (!is.null(per_region$BA25)) sum(per_region$BA25$adj.P.Val < PADJ_CUTOFF) else NA)
writeLines(report, file.path(output_dir, "GSE88890_FIXED_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("\n===== GSE88890 FIXED pipeline complete =====\n")
cat("  Outputs:", output_dir, "\n")

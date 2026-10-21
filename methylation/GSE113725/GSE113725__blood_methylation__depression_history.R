# ==============================================================================
# GSE113725__blood_methylation__depression_history.R (archived analysis copy)
# ------------------------------------------------------------------------------
# Dataset : GSE113725 - whole blood 450K DNA methylation, 194 samples
#           100 with self-reported history of depression / 100 without
#           (self-reported, NOT clinical MDD diagnosis - keep this caveat)
# Platform: GPL13534 Illumina HumanMethylation450
# Original: Crawford et al. 2018, Hum Mol Genet 27:2840 (PMID 29790996)
#           -> 6 Sidak-significant DMRs, top = LTB4R2 exon1 (P=1.27e-14)
# ------------------------------------------------------------------------------
# Fixes vs old pipeline (see All_Pipelines_Evaluation_Report.md):
#   [P0] Dataset was REJECTED as "insufficient power", but the real cause was
#        a WRONG groupid fold. Official GEO groupid coding (from GSM pages):
#           1 = depression history + inflammatory disorder history
#           2 = depression history only
#           3 = inflammatory disorder history only
#           4 = healthy controls
#        CORRECT depression contrast: {1,2} vs {3,4}
#        (actual counts in this deposit: 98 depressed / 96 control = 194 total;
#         GEO summary says 100/100 but only 194 samples are present)
#   [P1] Now uses age + gender as covariates (available in characteristics).
#   [P1] Explicit LTB4R2 (chr14 exon1) reproduction check.
#   [P2] Single output dir: <project>/methylation/GSE113725/results_fixed/
# Caveat : phenotype is SELF-REPORTED history; do not merge with clinical MDD
#          cohorts without stating this.
# Thresholds (user policy): adj.P < 0.05; no beta-diff cutoff.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  # NOTE: IlluminaHumanMethylation450kanno.ilmn12.hg19 is loaded INSIDE
  # MODULE 2 with requireNamespace() so the script still runs when the
  # annotation package is missing (chr/SNP filtering is skipped then).
})

# -------------------------------------------------------------------------------
# SETTINGS
# -------------------------------------------------------------------------------
project_root <- "/Users/tan/Documents/MDD-MultiOmics/transcriptomics"
output_dir   <- file.path(project_root, "methylation/GSE113725/results_fixed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

series_file <- "/Users/tan/Developer/projects/r_project/Study/GSE113725/GSE113725_series_matrix.txt"
manifest_file <- "/Users/tan/Developer/projects/r_project/Study/GSE88890/GSE88890_RAW/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
PADJ_CUTOFF <- 0.05

# ==============================================================================
# MODULE 1: LOAD METADATA & DATA (groupid official coding - CORE FIX)
# ==============================================================================
cat("\n===== MODULE 1: LOAD METADATA & DATA =====\n")

# --- 1.1 Parse header (age / groupid / gender) --------------------------------
parse_sm <- function(f, n_meta_lines = 200) {
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
# NOTE: `lst` is the WHOLE row vector (194 values); `lst[[1]]` is only the
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

gsm_ids <- sm[["Sample_geo_accession"]][[1]]
age     <- suppressWarnings(as.numeric(extract_field(sm[["Sample_characteristics_ch1"]], "age")))
groupid <- suppressWarnings(as.integer(extract_field(sm[["Sample_characteristics_ch1"]], "groupid")))
gender  <- extract_field(sm[["Sample_characteristics_ch1"]], "gender")

cat("  groupid raw counts:\n"); print(table(groupid, useNA = "ifany"))

# --- 1.2 FIX [P0]: official groupid coding -------------------------------------
#   1 = depression + inflammatory disorder; 2 = depression only
#   3 = inflammatory disorder only;       4 = healthy control
Group <- ifelse(groupid %in% c(1, 2), "depressed",
         ifelse(groupid %in% c(3, 4), "control", NA))
pheno <- data.frame(Sample_ID = gsm_ids, Group = Group,
                    groupid = groupid, Age = age, Gender = as.factor(gender),
                    stringsAsFactors = FALSE)
cat("  Groupid distribution (official coding):\n"); print(table(pheno$groupid, pheno$Group, useNA = "ifany"))
cat("  Depression contrast: depressed", sum(pheno$Group == "depressed"),
    "vs control", sum(pheno$Group == "control"), "\n")

# --- 1.3 Read data table (series matrix table block; 450K probes x 194) --------
cat("  Reading beta table (this is a 1.5GB file, takes ~1-2 min)...\n")
tab <- fread(series_file, skip = "!series_matrix_table_begin", header = TRUE,
             sep = "\t", na.strings = "NA", showProgress = FALSE)
# strip the begin marker row
tab <- tab[!grepl("^!series_matrix_table_end", tab[[1]]), ]
probe_id <- tab[[1]]
beta_mat <- as.matrix(tab[, -1, with = FALSE])
rownames(beta_mat) <- probe_id
colnames(beta_mat) <- gsm_ids
cat(sprintf("  Beta matrix: %d probes x %d samples\n", nrow(beta_mat), ncol(beta_mat)))
rm(tab); gc()

# Align metadata
pheno <- pheno[match(colnames(beta_mat), pheno$Sample_ID), ]
stopifnot(!any(is.na(pheno$Sample_ID)))
pheno <- pheno[!is.na(pheno$Group), ]
beta_mat <- beta_mat[, pheno$Sample_ID]
cat("  Samples with valid Group:", nrow(pheno), "\n")

# ==============================================================================
# MODULE 2: QC & PROBE FILTERING
# ==============================================================================
cat("\n===== MODULE 2: QC & PROBE FILTERING =====\n")

# --- 2.1 Load 450K annotation (self-contained, fallback-safe) ------------------
# NOTE: `getAnnotation()` belongs to minfi. The annotation package exports an
# object with the SAME name, so as.data.frame(<pkg object>) works directly.
# If the package is not installed, chr/SNP filtering is SKIPPED (NA-rate filter
# always runs) and a clear warning is printed. Install once with:
#   BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  ann <- as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann$Probe <- rownames(ann)
  cat(sprintf("  Annotation loaded: %d probes\n", nrow(ann)))
} else if (file.exists(manifest_file)) {
  raw_ann <- fread(manifest_file, skip = 7, fill = TRUE, showProgress = FALSE,
                   select = c("Name", "CHR", "MAPINFO", "Probe_SNPs", "Probe_SNPs_10",
                              "UCSC_RefGene_Name", "UCSC_RefGene_Group",
                              "Relation_to_UCSC_CpG_Island"))
  raw_ann <- raw_ann[grepl("^(cg|ch|rs)", Name)]
  ann <- data.frame(
    Probe = raw_ann$Name,
    chr = ifelse(grepl("^chr", raw_ann$CHR), raw_ann$CHR, paste0("chr", raw_ann$CHR)),
    pos = as.numeric(raw_ann$MAPINFO),
    SNP_ID = paste(ifelse(is.na(raw_ann$Probe_SNPs), "", raw_ann$Probe_SNPs),
                   ifelse(is.na(raw_ann$Probe_SNPs_10), "", raw_ann$Probe_SNPs_10)),
    UCSC_RefGene_Name = raw_ann$UCSC_RefGene_Name,
    UCSC_RefGene_Group = raw_ann$UCSC_RefGene_Group,
    Relation_to_Island = raw_ann$Relation_to_UCSC_CpG_Island,
    stringsAsFactors = FALSE)
  rm(raw_ann)
  cat(sprintf("  Annotation loaded from local GPL13534 manifest: %d probes\n", nrow(ann)))
} else {
  ann <- NULL
  stop("450K annotation unavailable: install the annotation package or restore the local GPL13534 manifest.")
}

# NA-rate filter (always runs, needs no annotation)
keep_probe <- rowMeans(is.na(beta_mat)) <= 0.05
# sex-chr + SNP probes (only when annotation is available)
if (!is.null(ann)) {
  sex_chr <- ann$Probe[ann$chr %in% c("chrX", "chrY")]
  snp_probes <- if ("SNP_ID" %in% colnames(ann)) ann$Probe[!is.na(ann$SNP_ID) & nzchar(trimws(ann$SNP_ID))] else character(0)
  keep_probe <- keep_probe & !(rownames(beta_mat) %in% c(sex_chr, snp_probes))
}
beta_qc <- beta_mat[keep_probe, ]
cat(sprintf("  Probes kept: %d\n", nrow(beta_qc)))

# ==============================================================================
# MODULE 3: DMP ANALYSIS  (depressed {1,2} vs control {3,4}, adjusted)
# ==============================================================================
cat("\n===== MODULE 3: DMP ANALYSIS =====\n")
pheno$Group <- factor(pheno$Group, levels = c("control", "depressed"))
design <- model.matrix(~ Group + Age + Gender, data = pheno)
colnames(design) <- make.names(colnames(design))
cat("  Design columns:", colnames(design), "\n")
if (qr(design)$rank != ncol(design)) stop("Methylation design matrix is rank deficient.")

delta_beta <- rowMeans(beta_qc[, pheno$Group == "depressed", drop = FALSE], na.rm = TRUE) -
  rowMeans(beta_qc[, pheno$Group == "control", drop = FALSE], na.rm = TRUE)
m_qc <- log2(pmin(pmax(beta_qc, 1e-6), 1 - 1e-6) /
               (1 - pmin(pmax(beta_qc, 1e-6), 1 - 1e-6)))
fit <- lmFit(m_qc, design)
fit <- eBayes(fit)
dmp_all <- topTable(fit, coef = "Groupdepressed", number = Inf, adjust.method = "BH")
dmp_all$ProbeID <- rownames(dmp_all)
dmp_all$DeltaBeta <- delta_beta[dmp_all$ProbeID]
dmp_sig <- dmp_all[dmp_all$adj.P.Val < PADJ_CUTOFF, ]

cat(sprintf("  DMPs adj.P<0.05 (depressed vs control): %d\n", nrow(dmp_sig)))
cat("  Up-methylated:", sum(dmp_sig$logFC > 0), "| Down-methylated:", sum(dmp_sig$logFC < 0), "\n")
write.csv(dmp_all, file.path(output_dir, "GSE113725_DMP_all.csv"), row.names = FALSE)
write.csv(dmp_sig, file.path(output_dir, "GSE113725_DMP_significant_adjP0.05.csv"), row.names = FALSE)

pvals <- dmp_all$P.Value[is.finite(dmp_all$P.Value)]
lambda <- median(qchisq(1 - pvals, 1)) / qchisq(0.5, 1)
diag_df <- data.frame(
  Probes_tested = length(pvals), Min_P = min(pvals),
  P_lt_0.05 = sum(pvals < 0.05), Expected_P_lt_0.05 = 0.05 * length(pvals),
  FDR_lt_0.05 = sum(dmp_all$adj.P.Val < 0.05, na.rm = TRUE),
  Lambda = lambda)
write.csv(diag_df, file.path(output_dir, "GSE113725_pvalue_diagnostics.csv"), row.names = FALSE)
cat(sprintf("  P-value diagnostic: lambda=%.3f; P<0.05=%d (expected %.0f)\n",
            lambda, sum(pvals < 0.05), 0.05 * length(pvals)))

plot_dmp <- dmp_all[is.finite(dmp_all$DeltaBeta) & is.finite(dmp_all$P.Value) &
                      dmp_all$P.Value > 0, ]
pdf(file.path(output_dir, "GSE113725_volcano_delta_beta.pdf"), width = 8, height = 6)
plot(plot_dmp$DeltaBeta, -log10(plot_dmp$P.Value), pch = 16, cex = 0.35,
     col = ifelse(plot_dmp$adj.P.Val < PADJ_CUTOFF, "#D7301F", "grey72"),
     main = "GSE113725: depression history vs control",
     sub = "M-value inference; x-axis shows delta-beta; colour = BH FDR",
     xlab = "Delta beta (depressed - control)", ylab = "-log10(P)")
abline(h = -log10(0.05), lty = 2, col = "#2C7FB8")
dev.off()

# ==============================================================================
# MODULE 4: ANNOTATION & LTB4R2 CHECK (reproduce the paper's top DMR)
# ==============================================================================
cat("\n===== MODULE 4: ANNOTATION & LTB4R2 CHECK =====\n")

# DMP gene/chr annotation (only when the annotation package is available)
if (nrow(dmp_sig) > 0 && !is.null(ann)) {
  dmp_sig$Chr       <- ann$chr[match(dmp_sig$ProbeID, ann$Probe)]
  dmp_sig$Pos       <- ann$pos[match(dmp_sig$ProbeID, ann$Probe)]
  dmp_sig$Gene      <- ann$UCSC_RefGene_Name[match(dmp_sig$ProbeID, ann$Probe)]
  dmp_sig$GeneGroup <- ann$UCSC_RefGene_Group[match(dmp_sig$ProbeID, ann$Probe)]
  write.csv(dmp_sig, file.path(output_dir, "GSE113725_DMP_annotated_adjP0.05.csv"), row.names = FALSE)
} else if (nrow(dmp_sig) > 0) {
  cat("  [WARN] Annotation package missing; gene annotation skipped.\n")
}

# LTB4R2 = chr14:24,789,000-24,791,000 (hg19 approx, exon1 target)
ltb_hits <- data.frame()     # initialize so MODULE 5 can always reference it
if (!is.null(ann)) {
  ltb_probes <- ann$Probe[ann$chr == "chr14" & ann$pos > 24780000 & ann$pos < 24800000]
  ltb_hits <- dmp_all[dmp_all$ProbeID %in% ltb_probes, ]
  cat(sprintf("  LTB4R2-region probes tested: %d\n", length(ltb_probes)))
  if (nrow(ltb_hits) > 0) {
    cat("  LTB4R2 region results (paper: top DMR P=1.27e-14, Sidak):\n")
    print(ltb_hits[, c("ProbeID", "logFC", "DeltaBeta", "P.Value", "adj.P.Val")])
    write.csv(ltb_hits, file.path(output_dir, "GSE113725_LTB4R2_check.csv"), row.names = FALSE)
  } else {
    cat("  No LTB4R2-region probes significant at probe level.\n")
    cat("  (Paper reports DMR-level Sidak significance; DMRcate/comb-p needed.)\n")
  }
} else {
  cat("  [WARN] Annotation package missing; LTB4R2 region check skipped.\n")
}

# ==============================================================================
# MODULE 5: REPORT
# ==============================================================================
report <- sprintf(
"GSE113725 FIXED pipeline report (2026-09)
============================================
Samples       : %d (depressed=%d [groupid 1+2], control=%d [groupid 3+4])
FIX [P0]      : official groupid coding restored - old script folded groups
                wrongly and the dataset was rejected as 'underpowered'.
                Actual: groupid 1=49, 2=49, 3=48, 4=48 (98/96 contrast)
Covariates    : Age + Gender
Probes        : %d (after NA/sex-chr/SNP filter)
DMPs adj.P<0.05: %d
P-value lambda: %.3f (values >1 suggest residual confounding/inflation)
LTB4R2 check  : %s
Caveat        : phenotype = SELF-REPORTED depression history (Crawford 2018).
                Do not merge with clinical MDD cohorts without stating this.
                Blood-cell proportions, array/slide effects and cross-reactive
                probes were not available/handled in this run; re-fit with
                these controls before treating the result as definitive.
",
  nrow(pheno), sum(pheno$Group == "depressed"), sum(pheno$Group == "control"),
  nrow(beta_qc), nrow(dmp_sig), lambda,
  if (nrow(ltb_hits) > 0) sprintf("%d probes tested; min P=%.2e; FDR<0.05=%d (DMR test required)",
                                  nrow(ltb_hits), min(ltb_hits$P.Value, na.rm = TRUE),
                                  sum(ltb_hits$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE)) else "no annotated probes in window")
writeLines(report, file.path(output_dir, "GSE113725_FIXED_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("\n===== GSE113725 FIXED pipeline complete =====\n")
cat("  Outputs:", output_dir, "\n")

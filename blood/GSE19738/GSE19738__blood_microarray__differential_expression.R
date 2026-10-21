# ==============================================================================
# blood_GSE19738_pipeline.R
# PURPOSE  : Add ONE new peripheral-blood MDD cohort (GSE19738) to the blood
#            transcriptome layer. GSE19738 is Agilent GPL6848 whole blood from
#            non-medicated MDD patients and controls, measured both at baseline
#            and after ex-vivo LPS stimulation.  Primary validation uses ONLY
#            the 33 MDD / 34 control basal samples.
# INPUT    : none (downloads from GEO). Set geo_dir first.
# OUTPUT   : unified_tables/unified_deg_table_GSE19738.tsv
#            results/deg_GSE19738_<contrast>.csv
# DEPENDS  : BiocManager: GEOquery, limma, oligo (or affy); data.table
# NOTE     : Sample numbers on GEO: ~33 MDD / 34 control (verify on the GEO
#            page before interpreting). This is a MICROARRAY dataset - no
#            DESeq2; limma + eBayes on the normalized matrix from GEO.
# USAGE    : Rscript GSE19738__blood_microarray__differential_expression.R [geo_dir]
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(data.table)
})

# ------------------------------------------------------------------------------
# MODULE 0: CONFIG
# ------------------------------------------------------------------------------
# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
dataset_dir <- file.path(OUTPUT_ROOT, "blood/GSE19738")
unified_dir <- file.path(dataset_dir, "unified_tables")
res_dir <- file.path(dataset_dir, "results")
dir.create(unified_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
# GEO cache: prefer a pre-downloaded copy in the author workspace; otherwise a
# local cache under the output tree (GEOquery downloads it there if absent).
geo_default <- file.path(WORK_ROOT, "analysis_framework/GSE19738/cache")
if (!dir.exists(geo_default)) geo_default <- file.path(dataset_dir, "cache")
geo_dir <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else geo_default

PADJ_CUTOFF <- 0.05

# Shared helpers are defined locally so this server script is standalone.
clean_genes <- function(x) {
  x <- trimws(sub("///.*$", "", as.character(x)))
  x <- toupper(x)
  x[!grepl("^[A-Z0-9][A-Z0-9.\\-]*$", x)] <- NA_character_
  x
}

write_unified <- function(dt, id) {
  dt <- dt[!is.na(gene) & is.finite(logFC) & is.finite(P.Value)]
  dt[, Direction := ifelse(logFC > 0, "UP", "DOWN")]
  out <- file.path(unified_dir, paste0("unified_deg_table_", id, ".tsv"))
  fwrite(dt[, .(gene, logFC, P.Value, adj.P.Val, Direction)], out, sep = "\t")
  cat(sprintf("  Unified table written: %s (%d genes)\n", out, nrow(dt)))
}

cat("\n===== MODULE 0: CONFIG =====\n")
cat("  GEO cache dir:", geo_dir, "\n\n")

# ------------------------------------------------------------------------------
# MODULE 1: DOWNLOAD + METADATA
# ------------------------------------------------------------------------------
cat("===== MODULE 1: DOWNLOAD + METADATA =====\n")
series_cache <- file.path(geo_dir, "GSE19738_series_matrix.txt.gz")
if (file.exists(series_cache)) {
  cat("  Using cached series matrix:", series_cache, "\n")
  gse <- getGEO(filename = series_cache, getGPL = FALSE)
} else {
  gse <- getGEO("GSE19738", destdir = geo_dir, getGPL = FALSE)[[1]]
}
pd <- pData(gse)

# Characteristics: this series uses "characteristics_ch1" with fields such as
#   diagnosis: Major Depressive Disorder / Control
#   (verify actual labels on the GEO page)
ch <- grep("characteristics", colnames(pd), value = TRUE, ignore.case = TRUE)
char_list <- lapply(ch, function(cn) as.character(pd[[cn]]))

extract_char <- function(meta_list, label) {
  for (lst in meta_list) {
    if (length(lst) == 0) next
    field <- strsplit(lst[1], ": ", fixed = TRUE)[[1]][1]
    if (identical(field, label)) return(sub("^[^:]*: ", "", lst))
  }
  rep(NA_character_, length(meta_list[[1]]))
}

diag <- tryCatch(extract_char(char_list, "diagnosis"),
                 error = function(e) rep(NA, nrow(pd)))
if (all(is.na(diag))) diag <- extract_char(char_list, "group")
if (all(is.na(diag))) diag <- extract_char(char_list, "disease")
if (all(is.na(diag)) && "title" %in% names(pd)) diag <- as.character(pd$title)
age <- suppressWarnings(as.numeric(extract_char(char_list, "age")))
sex <- extract_char(char_list, "gender")
condition <- if ("title" %in% names(pd)) {
  ifelse(grepl("basal", pd$title, ignore.case = TRUE), "basal",
         ifelse(grepl("LPS", pd$title, ignore.case = TRUE), "LPS", NA_character_))
} else rep(NA_character_, nrow(pd))

cat("  Diagnosis distribution:\n"); print(table(diag, useNA = "ifany"))

mdd <- ifelse(grepl("mdd|depress|case", diag, ignore.case = TRUE), 1,
       ifelse(grepl("ctrl|control|healthy|normal", diag, ignore.case = TRUE), 0, NA))
if (sum(mdd == 1, na.rm = TRUE) == 0 || sum(mdd == 0, na.rm = TRUE) == 0) {
  cat("  [FATAL] Cannot resolve MDD/control groups.\n")
  cat("  Inspect the GEO page for GSE19738 and edit the regex above.\n")
  quit(status = 1)
}
cat(sprintf("  MDD: %d, Control: %d\n", sum(mdd == 1, na.rm = TRUE),
            sum(mdd == 0, na.rm = TRUE)))

# ------------------------------------------------------------------------------
# MODULE 2: EXPRESSION MATRIX
# ------------------------------------------------------------------------------
cat("\n===== MODULE 2: EXPRESSION MATRIX =====\n")
# Microarray: use the series matrix normalized expression if present
expr_mat <- exprs(gse)          # probes x samples (series-matrix normalized)
cat(sprintf("  Series matrix expression: %d probes x %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))

# Better: GEO supplementary RMA-normalized matrix is usually available;
# if the series matrix is absent/empty, download supplement.
if (nrow(expr_mat) == 0) {
  cat("  [INFO] Empty series-matrix expression; downloading supplement...\n")
  getGEOSuppFiles("GSE19738", baseDir = geo_dir)
  cand <- list.files(file.path(geo_dir, "GSE19738"),
                     pattern = "txt|csv|gz|CEL", full.names = TRUE)
  cat("  Candidates:\n"); print(cand)
  # If raw CEL files are present, process with oligo/affy:
  cel <- cand[grepl("CEL", cand, ignore.case = TRUE)]
  if (length(cel) > 0) {
    cat("  [INFO] Raw CEL files found - running RMA with oligo...\n")
    if (requireNamespace("oligo", quietly = TRUE)) {
      library(oligo)
      raw <- read.celfiles(cel)
      expr_mat <- exprs(rma(raw))
      colnames(expr_mat) <- sub("\\.CEL$", "", colnames(expr_mat))
    } else {
      stop("oligo not installed; install via BiocManager::install('oligo')")
    }
  } else {
    stop("No usable expression file found in supplements.")
  }
}

# Align phenotype to expression columns
col_ids <- colnames(expr_mat)
mdd <- mdd[match(col_ids, rownames(pd))]
age <- age[match(col_ids, rownames(pd))]
sex <- sex[match(col_ids, rownames(pd))]
condition <- condition[match(col_ids, rownames(pd))]
keep <- !is.na(mdd) & condition == "basal" &
  colSums(is.na(expr_mat)) < nrow(expr_mat) * 0.5
expr_mat <- expr_mat[, keep]
mdd <- mdd[keep]; age <- age[keep]; sex <- sex[keep]
cat(sprintf("  After alignment: %d samples (MDD=%d, Ctrl=%d)\n",
            ncol(expr_mat), sum(mdd == 1), sum(mdd == 0)))
if (ncol(expr_mat) != 67 || sum(mdd == 1) != 33 || sum(mdd == 0) != 34)
  stop("Expected the 67 basal samples (33 MDD / 34 control); check title parsing and GEO input.")

# Agilent sample tables contain literal `null` values for probes that failed
# the study's detection filter.  Carrying those sparse rows into robust eBayes
# produces NA prior weights.  The publication reports 22,789 retained basal
# spots; use complete finite basal probes before probe-to-gene collapse.
finite_probe <- rowSums(is.finite(expr_mat)) == ncol(expr_mat)
cat(sprintf("  Complete finite basal probes: %d / %d\n",
            sum(finite_probe), length(finite_probe)))
expr_mat <- expr_mat[finite_probe, , drop = FALSE]
if (nrow(expr_mat) < 10000) stop("Too few complete basal probes; check GEO parsing.")

# ------------------------------------------------------------------------------
# MODULE 3: GPL6848 ANNOTATION + GENE-LEVEL LIMMA
# ------------------------------------------------------------------------------
cat("\n===== MODULE 3: DIFFERENTIAL EXPRESSION =====\n")

# Fetch the platform actually used by this series (GPL6848, not GPL570).
platform_ids <- unique(as.character(pd$platform_id))
platform_ids <- platform_ids[!is.na(platform_ids) & nzchar(platform_ids)]
if (length(platform_ids) != 1) stop("Expected one GEO platform, found: ", paste(platform_ids, collapse = ", "))
cat("  Fetching platform annotation:", platform_ids, "\n")
gpl_cache <- file.path(geo_dir, paste0(platform_ids, ".soft.gz"))
gpl <- if (file.exists(gpl_cache)) {
  cat("  Using cached platform file:", gpl_cache, "\n")
  getGEO(filename = gpl_cache)
} else {
  getGEO(platform_ids, destdir = geo_dir)
}
gpl_tab <- Table(gpl)
gene_col <- grep("gene.symbol|gene.?symbol|symbol", colnames(gpl_tab),
                 ignore.case = TRUE, value = TRUE)[1]
if (is.na(gene_col)) stop("No gene-symbol column found in ", platform_ids)
probe2gene <- setNames(clean_genes(gpl_tab[[gene_col]]), gpl_tab$ID)
gene_symbol <- probe2gene[rownames(expr_mat)]
valid_gene <- !is.na(gene_symbol) & nzchar(gene_symbol)
probe_mean <- rowMeans(expr_mat[valid_gene, , drop = FALSE], na.rm = TRUE)
ord <- order(gene_symbol[valid_gene], -probe_mean)
gene_expr <- expr_mat[valid_gene, , drop = FALSE][ord, , drop = FALSE]
selected_symbol <- gene_symbol[valid_gene][ord]
selected_probe <- rownames(gene_expr)
first_gene <- !duplicated(selected_symbol)
gene_expr <- gene_expr[first_gene, , drop = FALSE]
selected_probe <- selected_probe[first_gene]
rownames(gene_expr) <- selected_symbol[first_gene]
cat(sprintf("  Gene-level matrix: %d genes x %d basal samples\n",
            nrow(gene_expr), ncol(gene_expr)))

design_cols <- list(MDD = factor(mdd))
if (mean(is.na(age)) < 0.5) design_cols$Age <- age
if (mean(is.na(sex)) < 0.5 && length(unique(na.omit(sex))) > 1)
  design_cols$Sex <- factor(sex)
model_df <- as.data.frame(design_cols)
complete_model <- complete.cases(model_df)
if (sum(!complete_model) > 0)
  cat(sprintf("  Removing %d sample(s) with missing model covariates.\n",
              sum(!complete_model)))
gene_expr <- gene_expr[, complete_model, drop = FALSE]
selected_probe <- selected_probe
model_df <- model_df[complete_model, , drop = FALSE]
rownames(model_df) <- colnames(gene_expr)
saveRDS(gene_expr, file.path(res_dir, "GSE19738_gene_expression.rds"))
saveRDS(model_df, file.path(res_dir, "GSE19738_metadata.rds"))
design <- model.matrix(as.formula(paste("~", paste(names(design_cols),
                                                   collapse = " + "))),
                       data = model_df)
if (nrow(design) != ncol(gene_expr))
  stop("Missing covariates caused samples to be dropped; resolve metadata before DE.")
if (qr(design)$rank != ncol(design)) stop("GSE19738 design matrix is rank deficient.")
fit <- lmFit(gene_expr, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
res <- topTable(fit, coef = "MDD1", number = Inf, sort.by = "none")
out <- data.table(ProbeID = selected_probe[match(rownames(res), rownames(gene_expr))],
                  gene = rownames(res),
                  logFC = res$logFC,
                  AveExpr = res$AveExpr,
                  t = res$t,
                  P.Value = res$P.Value,
                  adj.P.Val = res$adj.P.Val)
out <- out[!is.na(gene) & nzchar(gene)]

fwrite(out, file.path(res_dir, "deg_GSE19738.csv"))
write_unified(out, "GSE19738")
cat(sprintf("  GSE19738: %d genes tested, adj.P<0.05: %d\n",
            nrow(out), sum(out$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE)))

# Full-distribution volcano: raw P on the y-axis, BH FDR only for colour.
plot_deg <- out[is.finite(logFC) & is.finite(P.Value) & P.Value > 0]
pdf(file.path(res_dir, "GSE19738_volcano.pdf"), width = 8, height = 6)
plot(plot_deg$logFC, -log10(plot_deg$P.Value), pch = 16, cex = 0.55,
     col = ifelse(plot_deg$adj.P.Val < PADJ_CUTOFF, "#D7301F", "grey72"),
     main = "GSE19738 basal whole blood: MDD vs control",
     sub = sprintf("BH FDR <0.05: %d; dashed line = nominal P=0.05",
                   sum(out$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE)),
     xlab = "Adjusted log2 expression difference (MDD - control)",
     ylab = "-log10(P)")
abline(h = -log10(0.05), lty = 2, col = "#2C7FB8")
dev.off()

# Sex-stratified sensitivity analyses (same platform in both sexes here).
sex_outputs <- list(overall = out)
if ("Sex" %in% names(model_df)) {
  for (sx in c("female", "male")) {
    sx_idx <- tolower(as.character(model_df$Sex)) == sx
    sx_df <- droplevels(model_df[sx_idx, , drop = FALSE])
    sx_terms <- c("MDD", if ("Age" %in% names(sx_df)) "Age")
    sx_design <- model.matrix(as.formula(paste("~", paste(sx_terms, collapse = " + "))),
                              data = sx_df)
    if (qr(sx_design)$rank != ncol(sx_design)) stop("Rank-deficient ", sx, " design.")
    sx_fit <- eBayes(lmFit(gene_expr[, sx_idx, drop = FALSE], sx_design),
                       trend = TRUE, robust = TRUE)
    sx_res <- topTable(sx_fit, coef = "MDD1", number = Inf, sort.by = "none")
    sx_out <- data.table(
      ProbeID = selected_probe[match(rownames(sx_res), rownames(gene_expr))],
      gene = rownames(sx_res), logFC = sx_res$logFC, AveExpr = sx_res$AveExpr,
      t = sx_res$t, P.Value = sx_res$P.Value, adj.P.Val = sx_res$adj.P.Val)
    sx_out <- sx_out[!is.na(gene) & nzchar(gene)]
    sex_outputs[[sx]] <- sx_out
    fwrite(sx_out, file.path(res_dir, paste0("deg_GSE19738_", sx, ".csv")))
    cat(sprintf("  %s sensitivity: n=%d (MDD=%d, control=%d), FDR DEGs=%d\n",
                sx, sum(sx_idx), sum(as.character(sx_df$MDD) == "1"),
                sum(as.character(sx_df$MDD) == "0"),
                sum(sx_out$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE)))

    sx_plot <- sx_out[is.finite(logFC) & is.finite(P.Value) & P.Value > 0]
    pdf(file.path(res_dir, paste0("GSE19738_volcano_", sx, ".pdf")),
        width = 8, height = 6)
    plot(sx_plot$logFC, -log10(sx_plot$P.Value), pch = 16, cex = 0.55,
         col = ifelse(sx_plot$adj.P.Val < PADJ_CUTOFF, "#D7301F", "grey72"),
         main = paste("GSE19738", sx, "MDD vs control"),
         sub = sprintf("BH FDR <0.05: %d; dashed line = nominal P=0.05",
                       sum(sx_out$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE)),
         xlab = "Adjusted log2 expression difference (MDD - control)",
         ylab = "-log10(P)")
    abline(h = -log10(0.05), lty = 2, col = "#2C7FB8")
    dev.off()
  }
}

# ------------------------------------------------------------------------------
# MODULE 4: PRE-SPECIFIED OXIDATIVE-PHOSPHORYLATION VALIDATION
# ------------------------------------------------------------------------------
cat("\n===== MODULE 4: PRE-SPECIFIED GO:0006119 VALIDATION =====\n")
targets <- list()
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
  for (label in names(sex_outputs)) {
    current <- sex_outputs[[label]]
    ranks <- setNames(current$t, current$gene)
    ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
    gsea <- gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 1, seed = 1, eps = 0,
                  verbose = FALSE)
    gsea_df <- as.data.frame(gsea)
    suffix <- if (label == "overall") "" else paste0("_", label)
    fwrite(gsea_df, file.path(res_dir, paste0("GSE19738_GSEA_GO_BP", suffix, ".tsv")),
           sep = "\t")
    target <- gsea_df[gsea_df$ID == "GO:0006119", , drop = FALSE]
    targets[[label]] <- target
    fwrite(target,
           file.path(res_dir, paste0("GSE19738_predefined_GO0006119", suffix, ".tsv")),
           sep = "\t")
    if (nrow(target)) {
      cat(sprintf("  [%s] GO:0006119 OXPHOS: NES=%.3f, P=%.3g, BH q=%.3g\n",
                  label, target$NES[1], target$pvalue[1], target$p.adjust[1]))
    } else cat("  [", label, "] GO:0006119 was not returned.\n", sep = "")
  }
} else {
  cat("  clusterProfiler/org.Hs.eg.db unavailable; pre-specified GSEA skipped.\n")
}

target_value <- function(label, field) {
  target <- targets[[label]]
  if (!is.null(target) && nrow(target)) target[[field]][1] else NA_real_
}
summary_df <- rbindlist(lapply(names(sex_outputs), function(label) {
  z <- sex_outputs[[label]]
  data.table(Analysis = label, Genes_tested = nrow(z),
             FDR_DEGs = sum(z$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE),
             OXPHOS_NES = target_value(label, "NES"),
             OXPHOS_P = target_value(label, "pvalue"),
             OXPHOS_FDR = target_value(label, "p.adjust"))
}))
fwrite(summary_df, file.path(res_dir, "GSE19738_analysis_summary.tsv"), sep = "\t")
report <- sprintf(
"GSE19738 fixed local pipeline report
====================================
Samples: 67 basal whole-blood samples (MDD=33, control=34); non-medicated cohort
Complete finite basal probes: %d; gene-level tests: %d
DEGs at BH FDR <0.05: %d
GO:0006119 oxidative phosphorylation: NES=%.3f; P=%.3g; BH q=%.3g
Interpretation: this dataset is an independent pathway-level validation cohort;
absence of single-gene FDR hits does not negate a coordinated rank shift.

Sex-stratified sensitivity summary:
%s
",
  sum(finite_probe), nrow(out), sum(out$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE),
  target_value("overall", "NES"), target_value("overall", "pvalue"),
  target_value("overall", "p.adjust"),
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"))
writeLines(report, file.path(res_dir, "GSE19738_FIXED_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo.txt"))

cat("\n===== blood_GSE19738_pipeline.R complete =====\n")

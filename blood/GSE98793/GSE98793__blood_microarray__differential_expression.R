# ==============================================================================
# GSE98793_pipeline_fixed.R  (FIXED copy, 2026-09)
# ------------------------------------------------------------------------------
# Dataset : GSE98793 - whole blood transcriptome, MDD (n=128, 64 with GAD)
#           vs Control (n=64), GPL570 Affymetrix HG-U133 Plus 2.0
# Original: Leday et al. 2018, Transl Psychiatry 8:118 (PMID 28688579)
# ------------------------------------------------------------------------------
# Fixes vs old pipeline (see All_Pipelines_Evaluation_Report.md):
#   [P0-1] Immune BH was applied to only 5 post-hoc selected cell types
#          (double-dipping). NOW: BH over ALL deconvolution features.
#   [P0-2] Method/report mismatch (manuscript: MCP-counter+Wilcoxon;
#          old script: quanTIseq+xCell+Welch). NOW: MCP-counter as primary
#          method (matching manuscript) + quanTIseq secondary, single output dir.
#   [P1-1] Immune scores compared raw (no covariates). NOW: limma model
#          ~ disease + age + gender + batch + anxiety.
#   [P1-2] 64/128 MDD comorbid with GAD not separated. NOW: add GAD subgroup
#          sensitivity analysis (MDD_noGAD vs Control).
#   [P1-3] "medication sensitivity" is impossible (GEO has NO medication
#          field). NOW: removed; noted in report.
#   [P2-1] eBayes(fit2,...) referenced undefined fit2. NOW: fixed to fit.
#   [P2-2] Results scattered over Desktop/Psychology_study/工程目录. NOW:
#          single output dir: <project>/blood/GSE98793/results_fixed/
# Thresholds (user policy): significance = adj.P/padj < 0.05; no logFC cutoff.
# Module comments are in simple English. Main pipeline steps are kept.
# ==============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(matrixStats)
})

# -------------------------------------------------------------------------------
# SETTINGS (paths, thresholds)
# -------------------------------------------------------------------------------
project_root <- "/Users/tan/Documents/MDD-MultiOmics/transcriptomics"
data_dir     <- "/Users/tan/Developer/projects/r_project/Study/GSE98793/GSE98793_analysis"
output_dir   <- file.path(project_root, "blood/GSE98793/results_fixed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

PADJ_CUTOFF <- 0.05          # user policy: adj.P / padj < 0.05
LOGFC_CUT   <- NULL          # no logFC threshold (MDD effects are small)
mcpcounter_genes_file <- file.path(project_root, "reference/MCPcounter/genes.txt")

# ==============================================================================
# MODULE 1: LOAD DATA & PREPROCESSING (keep main steps of 01_preprocessing_qc.R)
# ==============================================================================
cat("\n===== MODULE 1: LOAD & PREPROCESSING =====\n")

# --- 1.1 Load series matrix (local first, GEO download fallback) -------------
local_file <- file.path(data_dir, "GSE98793_series_matrix.txt")
gse <- NULL
if (file.exists(local_file)) {
  tryCatch({
    # GPL annotation is supplied by hgu133plus2.db; avoid an unnecessary
    # network request that can make an otherwise local run fail.
    gse_list <- getGEO(filename = local_file, getGPL = FALSE)
    gse <- if (is.list(gse_list)) gse_list[[1]] else gse_list
  }, error = function(e) cat("  local load failed:", e$message, "\n"))
}
if (is.null(gse)) {
  gse_list <- getGEO("GSE98793", GSEMatrix = TRUE, getGPL = TRUE)
  gse <- gse_list[[1]]
}
expr_matrix <- exprs(gse)
pheno_data  <- pData(gse)
cat(sprintf("  Expression: %d probes x %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))

# --- 1.2 Map probes to gene symbols, collapse to genes (keep main step) ------
probe_ids <- rownames(expr_matrix)
mapped <- AnnotationDbi::select(hgu133plus2.db,
                                keys = probe_ids, keytype = "PROBEID",
                                columns = c("SYMBOL", "ENTREZID"))
mapped <- mapped[!duplicated(mapped$PROBEID), ]
sym <- setNames(mapped$SYMBOL, mapped$PROBEID)

expr_df <- data.frame(ProbeID = probe_ids,
                      GeneSymbol = sym[probe_ids],
                      expr_matrix, check.names = FALSE, stringsAsFactors = FALSE)
expr_df <- expr_df[!is.na(expr_df$GeneSymbol), ]
expr_df$mean_expr <- rowMeans(expr_df[, -(1:2), drop = FALSE])
expr_df <- expr_df[order(expr_df$GeneSymbol, -expr_df$mean_expr), ]
expr_df <- expr_df[!duplicated(expr_df$GeneSymbol), ]   # keep max-mean probe per gene
# Keep only the original sample columns.  The previous code accidentally kept
# mean_expr as a 193rd "sample", which misaligned all phenotype rows.
sample_cols <- colnames(expr_matrix)
expr_clean <- as.matrix(expr_df[, sample_cols, drop = FALSE])
rownames(expr_clean) <- expr_df$GeneSymbol
stopifnot(ncol(expr_clean) == ncol(expr_matrix))
cat(sprintf("  Gene-level matrix: %d genes x %d samples\n", nrow(expr_clean), ncol(expr_clean)))

# --- 1.3 Parse phenotype: disease, age, gender, batch, anxiety ----------------
pheno_data$Sample <- rownames(pheno_data)
get_char <- function(pd, field, pattern) {
  col <- grep(field, colnames(pd), value = TRUE)[1]
  if (is.na(col)) return(rep(NA, nrow(pd)))
  vals <- pd[[col]]
  sub(paste0("^.*", pattern, "\\s*"), "", vals)
}
subject_group <- if ("subject group:ch1" %in% colnames(pheno_data)) {
  pheno_data[["subject group:ch1"]]
} else {
  pheno_data$characteristics_ch1
}
pheno_data$disease_group <- ifelse(grepl("CASE|MDD", subject_group, ignore.case = TRUE),
                                    "MDD", "Control")
pheno_data$age    <- as.numeric(if ("age:ch1" %in% names(pheno_data)) pheno_data[["age:ch1"]] else get_char(pheno_data, "age", "age:"))
pheno_data$gender <- if ("gender:ch1" %in% names(pheno_data)) pheno_data[["gender:ch1"]] else get_char(pheno_data, "gender", "gender:")
pheno_data$batch  <- if ("batch:ch1" %in% names(pheno_data)) pheno_data[["batch:ch1"]] else get_char(pheno_data, "batch", "batch:")
pheno_data$anxiety<- if ("anxiety:ch1" %in% names(pheno_data)) pheno_data[["anxiety:ch1"]] else get_char(pheno_data, "anxiety", "anxiety:")
# GAD comorbidity: parse from title or supplementary (0/1 for the 64 comorbid MDD)
pheno_data$GAD <- as.integer(tolower(trimws(pheno_data$anxiety)) %in%
                              c("1", "yes", "true", "gad", "anxiety"))

# Align samples (columns of expr_clean <-> rows of pheno)
pheno <- pheno_data[match(colnames(expr_clean), pheno_data$Sample), ]
if (anyNA(pheno$Sample)) stop("Expression columns do not align with phenotype rows.")
rownames(pheno) <- pheno$Sample
cat("  Group counts:\n"); print(table(pheno$disease_group))
cat("  GAD-comorbid within MDD:", sum(pheno$disease_group == "MDD" & pheno$GAD == 1), "\n")

# --- 1.4 QC plots (keep main step) ------------------------------------------
pdf(file.path(output_dir, "GSE98793_QC_report.pdf"), width = 10, height = 8)
is_logged <- median(expr_clean, na.rm = TRUE) < 10
expr_log <- if (!is_logged) log2(expr_clean + 1) else expr_clean
pca_res <- prcomp(t(expr_log), scale. = TRUE, center = TRUE)
pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
                     Group = pheno$disease_group)
print(ggplot(pca_df, aes(PC1, PC2, color = Group)) + geom_point(size = 2, alpha = 0.7) +
        theme_minimal() + ggtitle("PCA (GSE98793)"))
dev.off()
cat("  QC PDF saved.\n")

# Save preprocessed data (single source of truth)
saveRDS(expr_clean, file.path(output_dir, "GSE98793_expression_cleaned.rds"))
saveRDS(pheno,      file.path(output_dir, "GSE98793_metadata.rds"))
write.csv(pheno,    file.path(output_dir, "GSE98793_metadata.csv"), row.names = FALSE)

# ==============================================================================
# MODULE 2: DIFFERENTIAL EXPRESSION (limma, covariates; keep main steps)
# ==============================================================================
cat("\n===== MODULE 2: DIFFERENTIAL EXPRESSION (limma) =====\n")

# --- 2.1 Design: MDD vs Control adjusted for batch/gender/age (GAD flagged) ---
pheno$disease_group <- factor(pheno$disease_group, levels = c("Control", "MDD"))
design <- model.matrix(~ disease_group + batch + gender + age, data = pheno)
colnames(design) <- make.names(colnames(design))
if (nrow(design) != ncol(expr_clean)) stop("Design/expression sample count mismatch.")
if (qr(design)$rank != ncol(design)) stop("Primary design matrix is rank deficient.")
cat("  Design columns:", colnames(design), "\n")

# --- 2.2 limma fit (FIX: old code called eBayes(fit2,...) with undefined fit2) -
fit <- lmFit(expr_clean, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
deg_all <- topTable(fit, coef = "disease_groupMDD", number = Inf, adjust.method = "BH")
deg_all$GeneSymbol <- rownames(deg_all)

cat(sprintf("  DEGs with adj.P < 0.05: %d\n", sum(deg_all$adj.P.Val < PADJ_CUTOFF)))
deg_sig <- deg_all[deg_all$adj.P.Val < PADJ_CUTOFF, ]
write.csv(deg_all, file.path(output_dir, "GSE98793_all_DEG_results.csv"), row.names = FALSE)
write.csv(deg_sig, file.path(output_dir, "GSE98793_sig_DEGs_adjP0.05.csv"), row.names = FALSE)

# Volcano plot uses all finite tests.  Colour denotes BH FDR, while the dashed
# line is only the nominal P=0.05 reference.
plot_deg <- deg_all[is.finite(deg_all$logFC) & is.finite(deg_all$P.Value) &
                      deg_all$P.Value > 0, ]
plot_deg$Significance <- ifelse(plot_deg$adj.P.Val < PADJ_CUTOFF,
                                "FDR < 0.05", "FDR >= 0.05")
top_lab <- head(plot_deg[order(plot_deg$P.Value), ], 12)
p_volcano <- ggplot(plot_deg, aes(logFC, -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.55, size = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#2C7FB8") +
  geom_text(data = top_lab, aes(label = GeneSymbol), color = "black",
            check_overlap = TRUE, vjust = -0.6, size = 3) +
  scale_color_manual(values = c("FDR < 0.05" = "#D7301F", "FDR >= 0.05" = "grey72")) +
  labs(title = "GSE98793 whole blood: MDD vs control",
       subtitle = sprintf("%d genes at BH FDR < 0.05; dashed line = nominal P=0.05",
                          nrow(deg_sig)),
       x = "Adjusted log2 expression difference", y = expression(-log[10](italic(P))),
       color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(output_dir, "GSE98793_volcano.tiff"), p_volcano,
       width = 8, height = 6, dpi = 300, compression = "lzw")
ggsave(file.path(output_dir, "GSE98793_volcano.pdf"), p_volcano,
       width = 8, height = 6)

# --- 2.3 Sensitivity: MDD without GAD vs Control ------------------------------
pheno$subgroup <- ifelse(pheno$disease_group == "MDD" & pheno$GAD == 1,
                        "MDD_GAD", ifelse(pheno$disease_group == "MDD", "MDD_noGAD", "Control"))
pheno$subgroup <- factor(pheno$subgroup, levels = c("Control", "MDD_noGAD", "MDD_GAD"))
design2 <- model.matrix(~ subgroup + batch + gender + age, data = pheno)
colnames(design2) <- make.names(colnames(design2))
if (qr(design2)$rank != ncol(design2)) stop("GAD-sensitivity design matrix is rank deficient.")
fit2_base <- lmFit(expr_clean, design2)
fit2 <- eBayes(fit2_base, trend = TRUE, robust = TRUE)
deg_nogad <- topTable(fit2, coef = "subgroupMDD_noGAD", number = Inf, adjust.method = "BH")
write.csv(deg_nogad, file.path(output_dir, "GSE98793_DEG_MDDnoGAD_vs_Ctrl.csv"), row.names = TRUE)
cat(sprintf("  Sensitivity (MDD no-GAD vs Ctrl) DEGs adj.P<0.05: %d\n",
            sum(deg_nogad$adj.P.Val < PADJ_CUTOFF)))

deg_gad <- topTable(fit2, coef = "subgroupMDD_GAD", number = Inf, adjust.method = "BH")
write.csv(deg_gad, file.path(output_dir, "GSE98793_DEG_MDDGAD_vs_Ctrl.csv"), row.names = TRUE)
cat(sprintf("  Sensitivity (MDD GAD vs Ctrl) DEGs adj.P<0.05: %d\n",
            sum(deg_gad$adj.P.Val < PADJ_CUTOFF)))

contrast_gad <- makeContrasts(GAD_vs_noGAD = subgroupMDD_GAD - subgroupMDD_noGAD,
                              levels = design2)
fit_gad_diff <- contrasts.fit(fit2_base, contrast_gad)
fit_gad_diff <- eBayes(fit_gad_diff, trend = TRUE, robust = TRUE)
deg_gad_diff <- topTable(fit_gad_diff, coef = "GAD_vs_noGAD", number = Inf,
                         adjust.method = "BH")
write.csv(deg_gad_diff,
          file.path(output_dir, "GSE98793_DEG_MDDGAD_vs_MDDnoGAD.csv"), row.names = TRUE)
cat(sprintf("  Direct MDD-GAD vs MDD-noGAD DEGs adj.P<0.05: %d\n",
            sum(deg_gad_diff$adj.P.Val < PADJ_CUTOFF)))

# ==============================================================================
# MODULE 3: IMMUNE INFILTRATION  (CORE FIX: full-feature BH + covariate model)
# ==============================================================================
cat("\n===== MODULE 3: IMMUNE INFILTRATION (FIXED statistics) =====\n")

expr_mat <- as.matrix(expr_clean)

# --- 3.1 Primary method: MCP-counter (matches manuscript §3.5) ----------------
mcp_res <- NULL
if (requireNamespace("immunedeconv", quietly = TRUE)) {
  mcp_res <- tryCatch(immunedeconv::deconvolute(expr_mat, method = "mcp_counter"),
                      error = function(e) NULL)
}
if (is.null(mcp_res) && requireNamespace("MCPcounter", quietly = TRUE)) {
  cat("  immunedeconv unavailable; using MCPcounter directly.\n")
  mcp_res <- tryCatch({
    if (!file.exists(mcpcounter_genes_file)) stop("Local MCPcounter genes.txt is missing")
    marker_genes <- read.table(mcpcounter_genes_file, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE, check.names = FALSE,
                               colClasses = "character")
    res <- MCPcounter::MCPcounter.estimate(expr_mat, featuresType = "HUGO_symbols",
                                           genes = marker_genes)
    data.frame(cell_type = rownames(res), res, check.names = FALSE)
  }, error = function(e) {
    cat("  MCPcounter could not load its marker signatures:", e$message, "\n")
    NULL
  })
}
mcp_stats <- data.frame(cell_type = character(), padj = numeric())
mcp_adj <- data.frame(cell_type = character(), padj = numeric())

if (!is.null(mcp_res)) {
mcp_long <- mcp_res %>%
  pivot_longer(-cell_type, names_to = "Sample", values_to = "Abundance") %>%
  left_join(pheno[, c("Sample", "disease_group", "age", "gender", "batch", "anxiety", "GAD")],
            by = "Sample")

# --- 3.2 Secondary method: quanTIseq (blood-optimised) ------------------------
quant_res <- if (requireNamespace("immunedeconv", quietly = TRUE)) tryCatch(
  immunedeconv::deconvolute(expr_mat, method = "quantiseq", arrays = TRUE),
  error = function(e) NULL) else NULL
quant_long <- NULL
if (!is.null(quant_res)) {
  quant_long <- quant_res %>%
    pivot_longer(-cell_type, names_to = "Sample", values_to = "Abundance") %>%
    left_join(pheno[, c("Sample", "disease_group", "age", "gender", "batch", "anxiety", "GAD")],
              by = "Sample")
}

# --- 3.3 FIXED statistics: (a) full-feature BH on ALL cell types --------------
compare_cells_fullBH <- function(long_df, method_name) {
  results <- long_df %>%
    group_by(cell_type) %>%
    summarise(
      MDD_mean      = mean(Abundance[disease_group == "MDD"], na.rm = TRUE),
      Ctrl_mean     = mean(Abundance[disease_group == "Control"], na.rm = TRUE),
      fold_change   = MDD_mean / Ctrl_mean,
      logFC         = log2(MDD_mean / Ctrl_mean),
      p_value       = tryCatch(wilcox.test(Abundance[disease_group == "MDD"],
                                           Abundance[disease_group == "Control"])$p.value,
                               error = function(e) NA_real_),
      .groups = "drop")
  results$padj <- p.adjust(results$p_value, method = "BH")   # BH over ALL features
  results$method <- method_name
  results[order(results$padj), ]
}
mcp_stats  <- compare_cells_fullBH(mcp_long, "MCPcounter")
if (!is.null(quant_long)) quant_stats <- compare_cells_fullBH(quant_long, "quanTIseq")

cat("  MCP-counter full-feature BH results (top 6):\n")
print(head(mcp_stats, 6))
cat(sprintf("  MCP-counter cell types with padj < 0.05: %d\n",
            sum(mcp_stats$padj < PADJ_CUTOFF)))

write.csv(mcp_stats,  file.path(output_dir, "GSE98793_MCPcounter_fullBH_stats.csv"),  row.names = FALSE)
if (!is.null(quant_long)) write.csv(quant_stats, file.path(output_dir, "GSE98793_quanTIseq_fullBH_stats.csv"), row.names = FALSE)

# --- 3.4 Covariate-adjusted overall-MDD immune model --------------------------
# Do not include anxiety here: GAD occurs only inside the MDD group, so adding
# it changes disease_groupMDD into the MDD-noGAD-vs-control coefficient.  GAD is
# handled explicitly in the subgroup sensitivity analysis below.
adjust_immune_scores <- function(long_df, method_name) {
  scores <- long_df %>%
    pivot_wider(id_cols = Sample, names_from = cell_type, values_from = Abundance) %>%
    left_join(pheno[, c("Sample", "disease_group", "age", "gender", "batch")],
              by = "Sample") %>%
    filter(!is.na(disease_group))
  out <- list()
  for (ct in setdiff(colnames(scores), c("Sample", "disease_group", "age", "gender", "batch"))) {
    f <- as.formula(paste0("`", ct, "` ~ disease_group + age + gender + batch"))
    fit_i <- tryCatch(lm(f, data = scores), error = function(e) NULL)
    if (is.null(fit_i)) next
    s <- summary(fit_i)$coefficients
    if ("disease_groupMDD" %in% rownames(s)) {
      out[[ct]] <- data.frame(cell_type = ct, beta = s["disease_groupMDD", 1],
                              p_value = s["disease_groupMDD", 4], method = method_name)
    }
  }
  res <- do.call(rbind, out)
  res$padj <- p.adjust(res$p_value, method = "BH")
  res[order(res$padj), ]
}
mcp_adj <- adjust_immune_scores(mcp_long, "MCPcounter_overall_adj")
cat("  Covariate-adjusted overall-MDD immune model (top 6):\n"); print(head(mcp_adj, 6))
write.csv(mcp_adj, file.path(output_dir, "GSE98793_immune_adjCovariates.csv"), row.names = FALSE)

# --- 3.5 GAD sensitivity: MDD-noGAD vs Control for each immune feature --------
gad_sens <- function(long_df) {
  long_df$subgroup <- ifelse(long_df$disease_group == "MDD" & long_df$GAD == 1,
                             "MDD_GAD", ifelse(long_df$disease_group == "MDD", "MDD_noGAD", "Control"))
  long_df %>%
    filter(subgroup != "MDD_GAD") %>%
    group_by(cell_type) %>%
    summarise(p_value = tryCatch(wilcox.test(Abundance[subgroup == "MDD_noGAD"],
                                             Abundance[subgroup == "Control"])$p.value,
                                 error = function(e) NA_real_),
              .groups = "drop") %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    arrange(padj)
}
mcp_gad <- gad_sens(mcp_long)
cat("  GAD sensitivity (MDD-noGAD vs Ctrl), top 6:\n"); print(head(mcp_gad, 6))
write.csv(mcp_gad, file.path(output_dir, "GSE98793_immune_GADsensitivity.csv"), row.names = FALSE)

# --- 3.6 Boxplots for top features (keep main step, adapt thresholds) ---------
top_cells <- head(mcp_stats$cell_type, min(6, nrow(mcp_stats)))
plot_dat <- mcp_long %>% filter(cell_type %in% top_cells)
p <- ggplot(plot_dat, aes(x = cell_type, y = Abundance, fill = disease_group)) +
  geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.15, alpha = 0.4, size = 0.7) +
  scale_fill_manual(values = c(MDD = "#E41A1C", Control = "#377EB8")) +
  labs(title = "GSE98793 Immune Features (MDD vs Control, full-feature BH)",
       x = "Cell type", y = "MCP-counter score") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "GSE98793_immune_boxplot.pdf"), p, width = 9, height = 6)
} else {
  cat("  Neither immunedeconv nor MCPcounter is available; immune module skipped.\n")
}

# ==============================================================================
# MODULE 4: FUNCTIONAL ENRICHMENT (GSEA + GO/KEGG; keep main steps)
# ==============================================================================
cat("\n===== MODULE 4: FUNCTIONAL ENRICHMENT =====\n")
enrich_ok <- requireNamespace("clusterProfiler", quietly = TRUE)
if (enrich_ok) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  sig_genes <- rownames(deg_sig)
  if (length(sig_genes) >= 5) {
    ego <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = PADJ_CUTOFF)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0)
      write.csv(as.data.frame(ego), file.path(output_dir, "GSE98793_GO_BP_enrichment.csv"), row.names = FALSE)
    cat("  GO enrichment done.\n")
  } else {
    cat("  Too few significant genes for GO enrichment.\n")
  }
  # Pre-ranked GSEA (fgsea style via clusterProfiler gseGO)
  ranks <- setNames(deg_all$t, rownames(deg_all))
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  if (length(ranks) > 10) {
    gsea <- tryCatch(gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                           pAdjustMethod = "BH", pvalueCutoff = 1, seed = 1),
                     error = function(e) NULL)
    if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0)
      write.csv(as.data.frame(gsea), file.path(output_dir, "GSE98793_GSEA_BP.csv"), row.names = FALSE)
  }
} else {
  cat("  clusterProfiler not installed; enrichment skipped (optional).\n")
}

# ==============================================================================
# SUMMARY (single-source report)
# ==============================================================================
report <- sprintf(
"GSE98793 FIXED pipeline report (2026-09)
==========================================
Input   : %d genes x %d samples (MDD=%d [%d with GAD], Control=%d)
DEGs (adj.P<0.05)          : %d  (no logFC threshold, user policy)
DEGs MDD-noGAD vs Ctrl     : %d  (GAD sensitivity)
DEGs MDD-GAD vs Ctrl       : %d
DEGs MDD-GAD vs MDD-noGAD  : %d
MCP-counter cells padj<0.05: %d  (full-feature BH - NOT double-dipping)
  top cell types: %s
Covariate-adjusted immune  : %s
medication sensitivity     : NOT POSSIBLE - GEO GSE98793 has no medication field
                             (use drug-free dataset GSE19738 instead)
",
  nrow(expr_clean), ncol(expr_clean),
  sum(pheno$disease_group == "MDD"), sum(pheno$disease_group == "MDD" & pheno$GAD == 1),
  sum(pheno$disease_group == "Control"),
  sum(deg_all$adj.P.Val < PADJ_CUTOFF), sum(deg_nogad$adj.P.Val < PADJ_CUTOFF),
  sum(deg_gad$adj.P.Val < PADJ_CUTOFF), sum(deg_gad_diff$adj.P.Val < PADJ_CUTOFF),
  sum(mcp_stats$padj < PADJ_CUTOFF, na.rm = TRUE),
  if (nrow(mcp_stats)) paste(head(mcp_stats$cell_type, 5), collapse = ", ") else "not run",
  paste(if (exists("mcp_adj") && nrow(mcp_adj) > 0) head(mcp_adj$cell_type, 3) else "n/a", collapse = ", "))
writeLines(report, file.path(output_dir, "GSE98793_FIXED_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("\n===== GSE98793 FIXED pipeline complete =====\n")
cat("  Outputs:", output_dir, "\n")

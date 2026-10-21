# =============================================================================
# GSE54563 - FIXED pipeline copy
# Brain transcriptomics (MDD vs Control), Illumina HumanHT-12 V3.0 (GPL6947)
# Tissue: Anterior Cingulate Cortex (ACC); design: 25 MDD-Control matched pairs
#
# Reference (correct citation, GEO official): Chang LC et al.
# PLoS One. 2014;9(3):e90980. PMID: 24608543
#
# Main modules (same structure as original scripts):
#   M0 setup -> M1 phenotype -> M2 loading -> M3 NA removal -> M4 log2 check
#   -> M5 expression filter -> M6 sample QC -> M7 paired limma DEA
#   -> M8 DEG selection -> M9 annotation -> M10 dedup -> M11 paired Cohen's dz
#   -> M12 enrichment -> M13 save results -> M14 MOFA export
#
# Fixes vs original:
#   1. Filtering : unified rule - "log2 expression > 3 in >=10% of samples"
#                  (old code had a bug: used '> 8' on possibly-log2 data;
#                   two scripts used different thresholds)
#   2. DEG calls : adjusted P (BH, adj.P.Val < 0.05); no logFC cutoff
#   3. Annotation : illuminaHumanv3.db (matches GPL6947);
#                  old post-analysis script wrongly used illuminaHumanv4.db
#   4. Enrichment : uses adj.P < 0.05 genes (was raw P < 0.05)
#   5. Citation : corrected to Chang et al. 2014
# =============================================================================

# -----------------------------------------------------------------------------
# MODULE 0: SETUP
# -----------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(limma)
library(matrixStats)
library(ggplot2)
library(AnnotationDbi)
library(illuminaHumanv3.db)   # FIX: v3 db matches GPL6947 (HumanHT-12 V3.0)

# Data file candidates (edit paths as needed)
data_candidates <- c(
  "GSE54563_series_matrix.txt",
  "data/GSE54563_series_matrix.txt",
  "/Users/tan/Developer/projects/r_project/Study/GSE54563/GSE54563_series_matrix.txt"
)
find_file <- function(candidates) {
  for (p in candidates) if (file.exists(p)) return(p)
  stop("Data file not found. Tried:\n  ", paste(candidates, collapse = "\n  "))
}
series_file <- find_file(data_candidates)

project_root <- "/Users/tan/Documents/MDD-MultiOmics/transcriptomics"
output_dir <- file.path(project_root, "brain/GSE54563/results_fixed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# FIX: enrichment now uses adjusted P < 0.05 genes (was raw P < 0.05)
P_ADJ_CUTOFF <- 0.05

# Publication-style plot saver (TIFF + PDF)
save_sci_plot <- function(p, f, w = 8, h = 6) {
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  ggsave(filename = f, plot = p, device = "tiff", width = w, height = h,
         dpi = 300, compression = "lzw")
  ggsave(filename = gsub("\\.tiff", ".pdf", f), plot = p, device = "pdf",
         width = w, height = h)
}

cat("=== GSE54563 FIXED pipeline started ===\n")

# -----------------------------------------------------------------------------
# MODULE 1: PHENOTYPE PARSING (lines 29-30 of series matrix)
# -----------------------------------------------------------------------------
pheno_raw <- fread(series_file, skip = 28, sep = "\t", header = FALSE,
                   nrows = 2, data.table = FALSE, fill = TRUE)

sample_titles <- gsub('"', '', as.character(pheno_raw[1, -1]))
sample_geo_ids <- gsub('"', '', as.character(pheno_raw[2, -1]))

disease_group <- ifelse(grepl("^MDD matched", sample_titles, ignore.case = TRUE),
                        "MDD", "Control")
pair_number <- as.numeric(gsub(".*#(\\d+)$", "\\1", sample_titles))

pheno_54563 <- data.frame(
  sample_id = sample_geo_ids,
  title = sample_titles,
  disease = disease_group,
  pair = pair_number,
  stringsAsFactors = FALSE
)
cat(sprintf("Metadata: %d samples (%d MDD, %d Control, %d pairs)\n",
            nrow(pheno_54563),
            sum(pheno_54563$disease == "MDD"),
            sum(pheno_54563$disease == "Control"),
            length(unique(pheno_54563$pair))))

# -----------------------------------------------------------------------------
# MODULE 2: EXPRESSION MATRIX LOADING (line 67 onwards)
# -----------------------------------------------------------------------------
expr_raw <- fread(series_file, skip = 66, sep = "\t", header = TRUE,
                  data.table = FALSE, fill = TRUE)
probe_ids <- expr_raw[, 1]
expr_matrix <- as.matrix(expr_raw[, -1])
rownames(expr_matrix) <- probe_ids
storage.mode(expr_matrix) <- "numeric"

expr_matrix <- expr_matrix[, match(pheno_54563$sample_id, colnames(expr_matrix))]
stopifnot(!anyNA(colnames(expr_matrix)))
stopifnot(identical(colnames(expr_matrix), pheno_54563$sample_id))
cat(sprintf("Expression matrix: %d probes x %d samples\n",
            nrow(expr_matrix), ncol(expr_matrix)))

# -----------------------------------------------------------------------------
# MODULE 3: NA REMOVAL
# -----------------------------------------------------------------------------
na_rows <- rowSums(is.na(expr_matrix)) > 0
expr_clean <- expr_matrix[!na_rows, ]
cat(sprintf("Probes with NA removed: %d, remaining: %d\n",
            sum(na_rows), nrow(expr_clean)))

# -----------------------------------------------------------------------------
# MODULE 4: LOG2 CHECK (auto-detect)
# -----------------------------------------------------------------------------
max_val <- max(expr_clean, na.rm = TRUE)
if (max_val > 50) {
  cat(sprintf("Raw scale detected (max=%.2f), applying log2(x+1)...\n", max_val))
  expr_clean <- log2(expr_clean + 1)
} else {
  cat("Data already on log2 scale, skipped.\n")
}

# -----------------------------------------------------------------------------
# MODULE 5: LOW-EXPRESSION FILTER (unified rule)
# -----------------------------------------------------------------------------
# FIX: use the same rule as the MOFA script:
#      keep probes with log2 expression > 3 (raw ~8) in >= 10% of samples.
#      The old part1 script used 'rowSums(expr_clean > 8)' on possibly-log2
#      data, which was an inconsistent / too strict threshold.
keep <- rowSums(expr_clean > 3) >= 0.1 * ncol(expr_clean)
expr_filtered <- expr_clean[keep, ]
cat(sprintf("Probes kept after filter: %d (%.1f%%)\n",
            sum(keep), 100 * sum(keep) / nrow(expr_clean)))

# -----------------------------------------------------------------------------
# MODULE 6: SAMPLE CORRELATION QC
# -----------------------------------------------------------------------------
cor_mat <- cor(expr_filtered, use = "complete.obs")
mdd_idx <- which(pheno_54563$disease == "MDD")
ctrl_idx <- which(pheno_54563$disease == "Control")
if (length(mdd_idx) > 1) {
  cat(sprintf("MDD within-group correlation: %.4f\n",
              mean(cor_mat[mdd_idx, mdd_idx][lower.tri(cor_mat[mdd_idx, mdd_idx])])))
}
if (length(ctrl_idx) > 1) {
  cat(sprintf("Control within-group correlation: %.4f\n",
              mean(cor_mat[ctrl_idx, ctrl_idx][lower.tri(cor_mat[ctrl_idx, ctrl_idx])])))
}

# -----------------------------------------------------------------------------
# MODULE 7: PAIRED LIMMA DEA (design: group + pair)
# -----------------------------------------------------------------------------
group <- factor(pheno_54563$disease, levels = c("Control", "MDD"))
pair <- factor(pheno_54563$pair)

design <- model.matrix(~ group + pair)  # pair factor absorbs matching variation
colnames(design)[2] <- "MDD_vs_Control"
if (qr(design)$rank != ncol(design)) stop("Paired design matrix is rank deficient.")

fit <- lmFit(expr_filtered, design)
contrast_matrix <- makeContrasts(MDD_vs_Control = `MDD_vs_Control`, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# -----------------------------------------------------------------------------
# MODULE 8: DEG SELECTION (adjusted P < 0.05, no logFC cutoff)
# -----------------------------------------------------------------------------
deg <- topTable(fit2, coef = "MDD_vs_Control", number = Inf, adjust.method = "BH")
deg$probe_id <- rownames(deg)
cat(sprintf("Total probes tested: %d\n", nrow(deg)))

# -----------------------------------------------------------------------------
# MODULE 9: ANNOTATION (illuminaHumanv3.db - FIXED from v4)
# -----------------------------------------------------------------------------
probe_means <- rowMeans(expr_filtered, na.rm = TRUE)
symbols <- mapIds(illuminaHumanv3.db, keys = rownames(deg), column = "SYMBOL",
                  keytype = "PROBEID", multiVals = "first")
deg$GeneSymbol <- symbols
deg <- deg[!is.na(deg$GeneSymbol) & deg$GeneSymbol != "", ]
deg$AvgExpr <- probe_means[rownames(deg)]

# -----------------------------------------------------------------------------
# MODULE 10: GENE DEDUP (keep probe with highest mean expression)
# -----------------------------------------------------------------------------
deg <- deg[order(deg$AvgExpr, decreasing = TRUE), ]
deg <- deg[!duplicated(deg$GeneSymbol), ]
cat(sprintf("Unique genes after dedup: %d\n", nrow(deg)))

# Significant DEG set (adj.P < 0.05, no logFC filter - MDD effect is small)
deg_sig <- deg[deg$adj.P.Val < P_ADJ_CUTOFF, ]
cat(sprintf("Significant DEGs (adj.P < 0.05): %d\n", nrow(deg_sig)))

# -----------------------------------------------------------------------------
# MODULE 11: EFFECT SIZE (paired Cohen's dz)
# -----------------------------------------------------------------------------
# The differential-expression model is paired, so the effect size must also be
# paired.  Pair columns explicitly and standardise the within-pair differences.
pair_levels <- sort(unique(pheno_54563$pair))
ctrl_order <- match(pair_levels, pheno_54563$pair[ctrl_idx])
mdd_order  <- match(pair_levels, pheno_54563$pair[mdd_idx])
if (anyNA(ctrl_order) || anyNA(mdd_order)) stop("Incomplete MDD-control pair detected.")
ctrl_expr <- expr_filtered[, ctrl_idx[ctrl_order], drop = FALSE]
mdd_expr  <- expr_filtered[, mdd_idx[mdd_order], drop = FALSE]
paired_diff <- mdd_expr - ctrl_expr
paired_mean_diff <- rowMeans(paired_diff, na.rm = TRUE)
paired_sd_diff <- rowSds(paired_diff, na.rm = TRUE)
cohens_dz <- paired_mean_diff / paired_sd_diff
cohens_dz[!is.finite(cohens_dz)] <- NA_real_
deg$paired_mean_difference <- paired_mean_diff[rownames(deg)]
deg$paired_sd_difference <- paired_sd_diff[rownames(deg)]
deg$cohens_dz <- cohens_dz[rownames(deg)]
cat(sprintf("Mean |paired Cohen's dz|: %.4f\n", mean(abs(cohens_dz), na.rm = TRUE)))

# -----------------------------------------------------------------------------
# MODULE 12: FUNCTIONAL ENRICHMENT (GO BP, KEGG, Reactome)
# -----------------------------------------------------------------------------
# FIX: use adj.P < 0.05 genes (was raw P < 0.05) to reduce false positives.
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

genes <- deg_sig$GeneSymbol
cat(sprintf("Genes with adj.P < %.2f: %d\n", P_ADJ_CUTOFF, length(genes)))

if (length(genes) >= 10) {
  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat(sprintf("Mapped ENTREZ IDs: %d\n", nrow(entrez)))

  if (nrow(entrez) > 0) {
    # --- GO Biological Process ---
    cat("\n--- Running GO BP enrichment ---\n")
    go <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
    if (!is.null(go) && nrow(go) > 0) {
      write.csv(as.data.frame(go), file.path(output_dir, "GSE54563_GO_BP_Enrichment.csv"),
                row.names = FALSE)
      # FIX: restore enrichment plots (were missing in the fixed copy)
      save_sci_plot(dotplot(go, showCategory = 15),
                    file.path(output_dir, "GSE54563_GO_BP_Dotplot.tiff"))
      cat(sprintf("GO enriched terms: %d\n", nrow(go)))
    } else {
      cat("No significant GO terms.\n")
    }

    # --- KEGG ---
    cat("\n--- Running KEGG enrichment ---\n")
    options(timeout = 600)
    if (requireNamespace("KEGG.db", quietly = TRUE)) {
      kegg <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                         pvalueCutoff = 0.1, qvalueCutoff = 0.2,
                         use_internal_data = TRUE)
    } else {
      kegg <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                         pvalueCutoff = 0.1, qvalueCutoff = 0.2)
    }
    if (!is.null(kegg) && nrow(kegg) > 0) {
      kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      write.csv(as.data.frame(kegg_readable),
                file.path(output_dir, "GSE54563_KEGG_Enrichment.csv"), row.names = FALSE)
      save_sci_plot(barplot(kegg, showCategory = 15),
                    file.path(output_dir, "GSE54563_KEGG_Barplot.tiff"))
      cat(sprintf("KEGG enriched pathways: %d\n", nrow(kegg)))
    } else {
      cat("No significant KEGG pathways.\n")
    }

    # --- Reactome ---
    cat("\n--- Running Reactome enrichment ---\n")
    react <- enrichPathway(gene = entrez$ENTREZID, organism = "human",
                           pvalueCutoff = 0.05, readable = TRUE)
    if (!is.null(react) && nrow(react) > 0) {
      write.csv(as.data.frame(react),
                file.path(output_dir, "GSE54563_Reactome_Enrichment.csv"), row.names = FALSE)
      save_sci_plot(dotplot(react, showCategory = 15),
                    file.path(output_dir, "GSE54563_Reactome_Dotplot.tiff"))
      cat(sprintf("Reactome enriched pathways: %d\n", nrow(react)))
    } else {
      cat("No significant Reactome pathways.\n")
    }
  }
} else {
  cat("Too few genes (<10) for enrichment.\n")
}

# -----------------------------------------------------------------------------
# MODULE 12b: VOLCANO PLOT (all tested genes, not only significant)
# -----------------------------------------------------------------------------
# Plot finite results only.  Colour represents FDR; the dashed line is labelled
# as the nominal P=0.05 reference and is not a significance claim.
plot_deg <- deg[is.finite(deg$logFC) & is.finite(deg$P.Value) & deg$P.Value > 0, ]
plot_deg$Significance <- ifelse(plot_deg$adj.P.Val < P_ADJ_CUTOFF,
                                "FDR < 0.05", "FDR >= 0.05")
top_lab <- head(plot_deg[order(plot_deg$P.Value), ], 10)
p_volcano <- ggplot(plot_deg, aes(logFC, -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.65, size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#2C7FB8") +
  geom_text(data = top_lab, aes(label = GeneSymbol), color = "black",
            check_overlap = TRUE, vjust = -0.6, size = 3) +
  scale_color_manual(values = c("FDR < 0.05" = "#D7301F", "FDR >= 0.05" = "grey70")) +
  labs(title = "GSE54563 ACC: paired MDD vs control",
       subtitle = sprintf("%d genes at BH FDR < 0.05; dashed line = nominal P=0.05",
                          nrow(deg_sig)),
       x = "Mean paired difference (MDD - control, log2 scale)",
       y = expression(-log[10](italic(P))), color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
save_sci_plot(p_volcano, file.path(output_dir, "GSE54563_Volcano.tiff"), 8, 6)

# -----------------------------------------------------------------------------
# MODULE 13: SAVE RESULTS
# -----------------------------------------------------------------------------
write.csv(deg[, c("probe_id", "GeneSymbol", "logFC", "paired_mean_difference",
                  "paired_sd_difference", "cohens_dz", "P.Value", "adj.P.Val")],
          file.path(output_dir, "GSE54563_All_DEG_Results.csv"), row.names = FALSE)
write.csv(deg_sig, file.path(output_dir, "GSE54563_Significant_DEGs_adjP0.05.csv"),
          row.names = FALSE)

summary_df <- data.frame(
  Metric = c("Total Samples", "MDD", "Control", "Pairs",
             "Probes Analyzed", "Unique Genes",
             "Sig DEG (adj.P<0.05)", "Mean |paired Cohen's dz|"),
  Value = c(nrow(pheno_54563),
            sum(pheno_54563$disease == "MDD"),
            sum(pheno_54563$disease == "Control"),
            length(unique(pheno_54563$pair)),
            nrow(expr_filtered),
            nrow(deg),
            nrow(deg_sig),
            sprintf("%.4f", mean(abs(cohens_dz), na.rm = TRUE)))
)
write.csv(summary_df, file.path(output_dir, "GSE54563_Analysis_Summary.csv"), row.names = FALSE)
cat("Results saved to:", output_dir, "\n")

# -----------------------------------------------------------------------------
# MODULE 14: MOFA EXPORT (gene-level matrix + metadata)
# -----------------------------------------------------------------------------
mofa_dir <- file.path(output_dir, "MOFA_inputs")
dir.create(mofa_dir, showWarnings = FALSE, recursive = TRUE)

# Rebuild gene-level matrix from the filtered probes (dedup by highest mean expr)
expr_mapped <- as.data.frame(expr_filtered)
expr_mapped$ProbeID <- rownames(expr_filtered)
expr_mapped$GeneSymbol <- mapIds(illuminaHumanv3.db, keys = rownames(expr_filtered),
                                 column = "SYMBOL", keytype = "PROBEID",
                                 multiVals = "first")
expr_mapped <- expr_mapped[!is.na(expr_mapped$GeneSymbol), ]
expr_mapped$mean_expr <- rowMeans(expr_mapped[, !(colnames(expr_mapped) %in%
                                                    c("ProbeID", "GeneSymbol", "mean_expr"))])
expr_mapped <- expr_mapped[order(expr_mapped$GeneSymbol, -expr_mapped$mean_expr), ]
expr_mapped <- expr_mapped[!duplicated(expr_mapped$GeneSymbol), ]
expr_mapped$mean_expr <- NULL

gene_vec <- expr_mapped$GeneSymbol
expr_final <- as.matrix(expr_mapped[, !(colnames(expr_mapped) %in%
                                          c("ProbeID", "GeneSymbol"))])
rownames(expr_final) <- gene_vec
cat(sprintf("MOFA gene-level matrix: %d genes x %d samples\n",
            nrow(expr_final), ncol(expr_final)))

saveRDS(expr_final, file.path(mofa_dir, "GSE54563_expression_cleaned.rds"))
saveRDS(pheno_54563, file.path(mofa_dir, "GSE54563_metadata.rds"))
write.csv(expr_final, file.path(mofa_dir, "GSE54563_expression_cleaned.csv"))

cat("\n=== GSE54563 FIXED pipeline complete ===\n")

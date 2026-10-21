# =============================================================================
# GSE80655 - FIXED pipeline copy (single integrated script)
# Brain RNA-seq (MDD vs Control), Illumina HiSeq 2000 (GPL11154)
# Regions: AnCg, DLPFC, nAcc; raw counts GRCh38.p13 (NCBI)
#
# Main modules (same structure as original 01_loader -> 04_quality_control):
#   M0 setup -> M1 load counts + ID map -> M2 phenotype -> M3 QC plots
#   -> M4 DESeq2 VST normalization -> M5 per-region DEA -> M6 Cohen's d
#   -> M7 enrichment (per region) -> M8 visualization -> M9 MOFA export
#
# Fixes vs original:
#   1. Removed dead duplicate code block in the region DEA loop
#   2. DEG calls : adjusted P (padj < 0.05); no logFC cutoff (MDD effect small)
#   3. Enrichment : run separately per brain region (was merged across regions),
#                   using padj < 0.05 genes
#   4. Phenotype : local series-matrix parsing first (no network dependency),
#                  GEOquery as fallback
#   5. Gene filter : at least 10 samples with >=1 count (standard DESeq2 pre-filter)
#   6. MA plot : average expression computed per region (was whole-dataset)
#   7. Batch note : design keeps ~ Diagnosis; add ~ Batch + Diagnosis if a
#                   batch column exists in the metadata
# =============================================================================

# -----------------------------------------------------------------------------
# MODULE 0: SETUP
# -----------------------------------------------------------------------------
library(DESeq2)
library(apeglm)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(RColorBrewer)

# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

# Data file candidates (edit paths as needed)
count_candidates <- c(
  "GSE80655_raw_counts_GRCh38.p13_NCBI.tsv",
  "data/GSE80655_raw_counts_GRCh38.p13_NCBI.tsv",
  file.path(DATA_ROOT, "Psychology_study/datasets/GSE80655_raw_counts_GRCh38.p13_NCBI.tsv")
)
series_candidates <- c(
  "GSE80655_series_matrix.txt",
  "data/GSE80655_series_matrix.txt",
  file.path(DATA_ROOT, "Psychology_study/datasets/GSE80655_series_matrix.txt")
)
metadata_candidates <- c(
  "data/GSE80655_metadata.csv",
  file.path(WORK_ROOT, "transcriptomics/brain/GSE80655/data/GSE80655_metadata.csv"),
  # Recovery path for a previous run, whose relative output directory
  # accidentally placed this file under GSE53987.
  file.path(WORK_ROOT, "transcriptomics/brain/GSE53987/scripts/results_fixed/MOFA_inputs/GSE80655_metadata.rds")
)
find_file <- function(candidates) {
  for (p in candidates) if (file.exists(p)) return(p)
  stop("Data file not found. Tried:\n  ", paste(candidates, collapse = "\n  "))
}
count_file <- find_file(count_candidates)
series_file <- if (any(file.exists(series_candidates))) find_file(series_candidates) else NULL

dataset_dir <- file.path(OUTPUT_ROOT, "brain/GSE80655")
output_dir <- file.path(dataset_dir, "results_fixed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# FIX: DEG significance threshold - adjusted P < 0.05, no logFC cutoff.
PADJ_CUTOFF <- 0.05

cat("=== GSE80655 FIXED pipeline started ===\n")

# -----------------------------------------------------------------------------
# MODULE 1: LOAD RAW COUNTS + ENTREZID -> SYMBOL MAP + DEDUP
# -----------------------------------------------------------------------------
GSE80655_data <- read.delim(file = count_file, header = TRUE, sep = "\t",
                            row.names = 1, check.names = FALSE)
cat(sprintf("Raw matrix: %d rows x %d columns\n", nrow(GSE80655_data), ncol(GSE80655_data)))

gene_ids <- rownames(GSE80655_data)
symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL",
                  keytype = "ENTREZID", multiVals = "first")
cat(sprintf("Mapped %d / %d genes to SYMBOL\n", sum(!is.na(symbols)), length(symbols)))

mapped_idx <- !is.na(symbols)
valid_symbols <- symbols[mapped_idx]
valid_matrix <- GSE80655_data[mapped_idx, ]

# Dedup: keep the gene with highest total count per symbol (same as original)
gene_sums <- rowSums(valid_matrix)
symbol_info <- data.frame(symbol = valid_symbols, sum_count = gene_sums,
                          index = seq_len(nrow(valid_matrix)))
symbol_info_sorted <- symbol_info[order(symbol_info$sum_count, decreasing = TRUE), ]
unique_info <- symbol_info_sorted[!duplicated(symbol_info_sorted$symbol), ]
final_matrix_unique <- valid_matrix[unique_info$index, ]
rownames(final_matrix_unique) <- unique_info$symbol
cat(sprintf("After dedup: %d unique gene symbols\n", nrow(final_matrix_unique)))

# -----------------------------------------------------------------------------
# MODULE 2: PHENOTYPE (local series matrix first, GEOquery fallback)
# -----------------------------------------------------------------------------
# FIX: parse the local series matrix (no network dependency).
#      Sample title format: X1834_AnCg_C_SL31501 -> [subject]_[region]_[dx]_[lib]
#      dx codes: C=Control, M=MDD, B=Bipolar, S=Schizophrenia
metadata_loaded <- FALSE

# Prefer a local, auditable metadata table when available.
metadata_file <- metadata_candidates[file.exists(metadata_candidates)][1]
if (!is.na(metadata_file)) {
  sample_metadata <- if (grepl("\\.rds$", metadata_file, ignore.case = TRUE)) {
    readRDS(metadata_file)
  } else {
    read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  }
  required_meta <- c("GSM_ID", "Sample_Title", "Diagnosis")
  metadata_loaded <- all(required_meta %in% colnames(sample_metadata))
  if (metadata_loaded) {
    sample_metadata <- sample_metadata[
      sample_metadata$Diagnosis %in% c("MDD", "Control"), , drop = FALSE]
    dir.create(file.path(dataset_dir, "data"), recursive = TRUE, showWarnings = FALSE)
    write.csv(sample_metadata,
              file.path(dataset_dir, "data/GSE80655_metadata.csv"), row.names = FALSE)
    cat(sprintf("Local metadata table loaded: %d MDD/control samples\n",
                nrow(sample_metadata)))
  }
}

if (!metadata_loaded && !is.null(series_file)) {
  tryCatch({
    lines <- readLines(series_file, n = 100)
    title_line <- lines[grep("^!Sample_title", lines)]
    gsm_line  <- lines[grep("^!Sample_geo_accession", lines)]
    if (length(title_line) > 0 && length(gsm_line) > 0) {
      titles <- gsub('"', '', strsplit(title_line[1], "\t")[[1]][-1])
      gsms   <- gsub('"', '', strsplit(gsm_line[1], "\t")[[1]][-1])

      parse_title <- function(t) {
        parts <- strsplit(t, "_")[[1]]
        dx_code <- parts[3]
        dx <- switch(dx_code, C = "Control", M = "MDD", B = "BD", S = "SZ", NA_character_)
        c(dx = dx)
      }
      diagnosis <- vapply(titles, function(t) unname(parse_title(t)["dx"]), character(1))

      sample_metadata <- data.frame(
        GSM_ID = gsms,
        Sample_Title = titles,
        Diagnosis = diagnosis,
        stringsAsFactors = FALSE
      )
      sample_metadata <- sample_metadata[!is.na(sample_metadata$Diagnosis) &
                                           sample_metadata$Diagnosis %in% c("MDD", "Control"), ]
      metadata_loaded <- nrow(sample_metadata) > 0
      if (metadata_loaded) cat(sprintf("Local metadata parsed: %d MDD/Control samples\n",
                                       nrow(sample_metadata)))
    }
  }, error = function(e) cat("Local metadata parse failed:", e$message, "\n"))
}

if (!metadata_loaded) {
  cat("Falling back to GEOquery (needs network)...\n")
  library(GEOquery)
  gse_data <- getGEO("GSE80655", GSEMatrix = TRUE)
  gse_eset <- gse_data[[1]]
  phenotype_df <- pData(gse_eset)
  phenotype_selected <- phenotype_df[, c("geo_accession", "brain region:ch1",
                                         "clinical diagnosis:ch1")]
  metadata_selected <- phenotype_selected[
    phenotype_selected$`clinical diagnosis:ch1` %in% c("Control", "Major Depression"), ]
  # FIX: also keep the sample title so Brain_Region can be parsed later.
  gsm_rownames <- rownames(metadata_selected)
  sample_metadata <- data.frame(
    GSM_ID = metadata_selected$geo_accession,
    Sample_Title = phenotype_df[gsm_rownames, "title"],
    Diagnosis = ifelse(metadata_selected$`clinical diagnosis:ch1` == "Control",
                       "Control", "MDD"),
    stringsAsFactors = FALSE
  )
  # add optional batch column if present
  if ("batch:ch1" %in% colnames(phenotype_df)) {
    sample_metadata$Batch <- phenotype_df[gsm_rownames, "batch:ch1"]
    cat("NOTE: batch:ch1 found -> Batch column added; DEA design will use ~ Batch + Diagnosis\n")
  }
  metadata_loaded <- TRUE
  dir.create(file.path(dataset_dir, "data"), recursive = TRUE, showWarnings = FALSE)
  write.csv(sample_metadata,
            file.path(dataset_dir, "data/GSE80655_metadata.csv"), row.names = FALSE)
}

cat("Diagnosis distribution:\n"); print(table(sample_metadata$Diagnosis))

# -----------------------------------------------------------------------------
# MODULE 2b: ALIGN COUNTS WITH METADATA, CONVERT TO INTEGER
# -----------------------------------------------------------------------------
selected_gsm_ids <- sample_metadata$GSM_ID
matching_gsm_ids <- selected_gsm_ids[selected_gsm_ids %in% colnames(final_matrix_unique)]
cat(sprintf("Samples matched to count matrix: %d / %d\n",
            length(matching_gsm_ids), length(selected_gsm_ids)))

expr_data_final <- as.matrix(final_matrix_unique[, matching_gsm_ids])
sample_metadata_final <- sample_metadata[match(colnames(expr_data_final),
                                               sample_metadata$GSM_ID), ]
if (anyNA(sample_metadata_final$GSM_ID)) stop("Count columns and metadata are misaligned.")
rownames(sample_metadata_final) <- sample_metadata_final$GSM_ID

expr_data_counts <- round(expr_data_final)
storage.mode(expr_data_counts) <- "integer"
cat(sprintf("Aligned integer counts: %d genes x %d samples\n",
            nrow(expr_data_counts), ncol(expr_data_counts)))

# Region info from sample titles (for per-region analysis)
parse_region <- function(t) {
  parts <- strsplit(t, "_")[[1]]
  r <- parts[2]
  ifelse(r %in% c("AnCg", "DLPFC", "nAcc"), r, NA_character_)
}
sample_metadata_final$Brain_Region <- vapply(sample_metadata_final$Sample_Title,
                                             parse_region, character(1))
regions <- unique(sample_metadata_final$Brain_Region)
regions <- regions[!is.na(regions)]
if (length(regions) == 0) {
  stop("Could not parse any brain region from sample titles. Check title format.")
}
regions_list <- list()
for (region in regions) {
  regions_list[[region]] <- sample_metadata_final$GSM_ID[
    sample_metadata_final$Brain_Region == region]
}
cat("Samples per region:\n"); print(table(sample_metadata_final$Brain_Region))

# -----------------------------------------------------------------------------
# MODULE 3: QC VISUALIZATION (sample correlation + PCA, pre-normalization)
# -----------------------------------------------------------------------------
save_sci_plot <- function(p, f, w = 8, h = 6) {
  ggsave(filename = file.path(output_dir, f), plot = p, device = "tiff",
         width = w, height = h, dpi = 300, compression = "lzw")
  ggsave(filename = gsub("\\.tiff", ".pdf", file.path(output_dir, f)),
         plot = p, device = "pdf", width = w, height = h)
}

log_expr <- log2(expr_data_counts + 1)

cor_matrix <- cor(log_expr, use = "complete.obs")
annotation_col <- data.frame(
  Diagnosis = sample_metadata_final$Diagnosis,
  Brain_Region = sample_metadata_final$Brain_Region,
  row.names = colnames(log_expr)
)
p_heatmap <- pheatmap(cor_matrix, annotation_col = annotation_col,
                      main = "Sample Correlation (GSE80655)", fontsize = 6, silent = TRUE)
tiff(file.path(output_dir, "QC_Sample_Correlation_Heatmap.tiff"),
     width = 10, height = 8, units = "in", res = 300, compression = "lzw")
grid::grid.draw(p_heatmap$gtable)
dev.off()

gene_var <- apply(log_expr, 1, var, na.rm = TRUE)
log_expr_var <- log_expr[gene_var > 0, ]
pca_res <- prcomp(t(log_expr_var), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$Diagnosis <- sample_metadata_final$Diagnosis
pca_df$Brain_Region <- sample_metadata_final$Brain_Region
var_explained <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis, shape = Brain_Region)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of GSE80655 Samples (pre-VST)",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("Control" = "#2E86AB", "MDD" = "#A23B72"))
save_sci_plot(p_pca, "QC_PCA.tiff")

# -----------------------------------------------------------------------------
# MODULE 4: DESEQ2 NORMALIZATION (size factors + VST, blind)
# -----------------------------------------------------------------------------
# FIX: standard pre-filter - gene detected (>=1 count) in at least 10 samples.
if ("Batch" %in% colnames(sample_metadata_final)) {
  design_formula <- ~ Batch + Diagnosis
} else {
  design_formula <- ~ Diagnosis
}
dds_all <- DESeqDataSetFromMatrix(countData = expr_data_counts,
                                  colData = sample_metadata_final,
                                  design = design_formula)
keep_genes <- rowSums(counts(dds_all) >= 1) >= 10
dds_all <- dds_all[keep_genes, ]
cat(sprintf("Genes kept after pre-filter (>=1 count in >=10 samples): %d\n", sum(keep_genes)))

dds_all <- estimateSizeFactors(dds_all)
vsd <- vst(dds_all, blind = TRUE)
vsd_matrix <- assay(vsd)
cat(sprintf("VST matrix: %d genes x %d samples\n", nrow(vsd_matrix), ncol(vsd_matrix)))

# Post-VST PCA (QC)
pca_data <- plotPCA(vsd, intgroup = c("Diagnosis", "Brain_Region"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 2)
p_pca_vst <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Diagnosis, shape = Brain_Region)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA (VST transformed) - GSE80655",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("Control" = "#2E86AB", "MDD" = "#A23B72"))
save_sci_plot(p_pca_vst, "QC_PCA_VST.tiff")

# -----------------------------------------------------------------------------
# MODULE 5: PER-REGION DIFFERENTIAL EXPRESSION (DESeq2 + apeglm)
# -----------------------------------------------------------------------------
# FIX: single clean loop (removed the dead duplicate block from the original).
# FIX: significance = padj < 0.05; no logFC cutoff.
# FIX: keep the FULL per-region result for volcano / MA plots (all tested genes).
deg_results_by_region <- list()
region_full_results <- list()
cohens_d_by_region <- list()

for (region in regions) {
  cat("\n>>> Processing region:", region, "\n")

  current_gsm <- regions_list[[region]]
  current_counts <- expr_data_counts[, current_gsm, drop = FALSE]
  current_metadata <- sample_metadata_final[match(current_gsm,
                                                  sample_metadata_final$GSM_ID), ]

  current_metadata$Diagnosis <- relevel(factor(current_metadata$Diagnosis), ref = "Control")

  if ("Batch" %in% colnames(current_metadata)) {
    dds <- DESeqDataSetFromMatrix(countData = current_counts,
                                  colData = current_metadata,
                                  design = ~ Batch + Diagnosis)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = current_counts,
                                  colData = current_metadata,
                                  design = ~ Diagnosis)
  }
  # Apply the same explicit pre-filter independently in each region.
  keep_region <- rowSums(counts(dds) >= 1) >= 10
  dds <- dds[keep_region, ]
  dds <- DESeq(dds, quiet = TRUE)

  # Use unshrunken Wald statistics for inference.  apeglm supplies a separate
  # effect-size estimate and requires a coefficient name, not a contrast.
  res_raw <- results(dds, contrast = c("Diagnosis", "MDD", "Control"))
  coef_name <- grep("Diagnosis.*MDD.*Control", resultsNames(dds), value = TRUE)[1]
  res_shrunk <- if (!is.na(coef_name)) tryCatch(
    lfcShrink(dds, coef = coef_name, type = "apeglm"),
    error = function(e) NULL) else NULL

  # FIX: padj < 0.05 (was pvalue < 0.01); no logFC cutoff
  # Full table for plotting (all tested genes)
  res_full <- as.data.frame(res_raw) %>%
    tibble::rownames_to_column("GeneSymbol") %>%
    dplyr::rename(logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj) %>%
    dplyr::mutate(Region = region)
  if (!is.null(res_shrunk)) {
    res_full$logFC_shrunk <- res_shrunk$log2FoldChange[
      match(res_full$GeneSymbol, rownames(res_shrunk))]
  } else {
    res_full$logFC_shrunk <- NA_real_
  }
  region_full_results[[region]] <- res_full
  write.csv(res_full,
            file.path(output_dir, paste0("GSE80655_all_DE_results_", region, ".csv")),
            row.names = FALSE)

  # Significant subset (for DEG CSV and enrichment)
  res_df <- res_full %>%
    dplyr::filter(adj.P.Val < PADJ_CUTOFF) %>%
    dplyr::select(GeneSymbol, logFC, P.Value, adj.P.Val, Region)

  cat(sprintf("  DEGs (padj < 0.05): %d\n", nrow(res_df)))
  deg_results_by_region[[region]] <- res_df

  # --- MODULE 6: EFFECT SIZE (Cohen's d, from VST) ---
  region_vst <- vsd_matrix[, colnames(vsd_matrix) %in% current_gsm, drop = FALSE]
  region_vst <- region_vst[, match(current_gsm, colnames(region_vst)), drop = FALSE]
  ctrl_idx <- current_metadata$Diagnosis == "Control"
  mdd_idx <- current_metadata$Diagnosis == "MDD"
  if (sum(ctrl_idx) > 1 & sum(mdd_idx) > 1) {
    ctrl_expr <- region_vst[, ctrl_idx, drop = FALSE]
    mdd_expr <- region_vst[, mdd_idx, drop = FALSE]
    ctrl_mean <- rowMeans(ctrl_expr, na.rm = TRUE)
    mdd_mean <- rowMeans(mdd_expr, na.rm = TRUE)
    ctrl_sd <- apply(ctrl_expr, 1, sd, na.rm = TRUE)
    mdd_sd <- apply(mdd_expr, 1, sd, na.rm = TRUE)
    cohens_d <- (mdd_mean - ctrl_mean) / sqrt((ctrl_sd^2 + mdd_sd^2) / 2)
    cohens_d_by_region[[region]] <- cohens_d
    cat(sprintf("  Mean |Cohen's d|: %.4f\n", mean(abs(cohens_d), na.rm = TRUE)))
  }
}

# -----------------------------------------------------------------------------
# MODULE 5b: MERGE DEG RESULTS
# -----------------------------------------------------------------------------
all_degs <- dplyr::bind_rows(deg_results_by_region)
cat(sprintf("\nTotal merged DEGs (padj < 0.05): %d\n", nrow(all_degs)))
# Always write a schema-valid table, including when there are zero hits.
write.csv(all_degs, file.path(output_dir, "GSE80655_MDD_vs_Control_BrainRegion_DEGs.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# MODULE 7: FUNCTIONAL ENRICHMENT - SEPARATELY PER REGION
# -----------------------------------------------------------------------------
# FIX: per-region enrichment (was merged across regions, which mixed signals).
enrich_region <- function(genes, region) {
  if (length(genes) < 10) {
    cat(sprintf("  %s: too few genes (<10) for enrichment, skipped.\n", region))
    return(invisible(NULL))
  }
  entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez_ids) == 0) {
    cat(sprintf("  %s: no ENTREZ mapping, skipped.\n", region))
    return(invisible(NULL))
  }
  cat(sprintf("  %s: running GO BP / KEGG / Reactome (%d mapped genes)...\n",
              region, nrow(entrez_ids)))

  go <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(go) && nrow(go) > 0) {
    write.csv(as.data.frame(go),
              file.path(output_dir, paste0("GSE80655_GO_BP_Enrichment_", region, ".csv")),
              row.names = FALSE)
    cat(sprintf("  %s: GO terms = %d\n", region, nrow(go)))
  }

  options(timeout = 600)
  if (requireNamespace("KEGG.db", quietly = TRUE)) {
    kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa",
                       pvalueCutoff = 0.1, qvalueCutoff = 0.2, use_internal_data = TRUE)
  } else {
    kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa",
                       pvalueCutoff = 0.1, qvalueCutoff = 0.2)
  }
  if (!is.null(kegg) && nrow(kegg) > 0) {
    kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(kegg_readable),
              file.path(output_dir, paste0("GSE80655_KEGG_Enrichment_", region, ".csv")),
              row.names = FALSE)
    cat(sprintf("  %s: KEGG pathways = %d\n", region, nrow(kegg)))
  }

  react <- enrichPathway(gene = entrez_ids$ENTREZID, organism = "human",
                         pvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(react) && nrow(react) > 0) {
    write.csv(as.data.frame(react),
              file.path(output_dir, paste0("GSE80655_Reactome_Enrichment_", region, ".csv")),
              row.names = FALSE)
    cat(sprintf("  %s: Reactome pathways = %d\n", region, nrow(react)))
  }
}

for (region in regions) {
  region_genes <- deg_results_by_region[[region]]$GeneSymbol
  cat("\n--- Enrichment for region:", region, "---\n")
  enrich_region(region_genes, region)
}

# -----------------------------------------------------------------------------
# MODULE 8: VISUALIZATION (volcano + MA, ALL tested genes per region)
# -----------------------------------------------------------------------------
# FIX: plot the full gene distribution per region (grey = non-significant,
#      red = adj.P < 0.05), instead of only the few significant DEGs.
for (region in regions) {
  full_res <- region_full_results[[region]]
  if (is.null(full_res) || nrow(full_res) == 0) next

  full_res$Significant <- ifelse(!is.na(full_res$adj.P.Val) &
                                   full_res$adj.P.Val < PADJ_CUTOFF,
                                 "Significant", "Not Significant")
  n_sig <- sum(full_res$Significant == "Significant")

  plot_res <- full_res[is.finite(full_res$logFC) &
                         is.finite(full_res$P.Value) & full_res$P.Value > 0, ]
  n_omitted <- nrow(full_res) - nrow(plot_res)
  top_lab <- head(plot_res[order(plot_res$P.Value), ], 8)

  p <- ggplot(plot_res, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = Significant), alpha = 0.6, size = 1.0) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#2C7FB8") +
    geom_text(data = top_lab, aes(label = GeneSymbol), color = "black",
              check_overlap = TRUE, vjust = -0.6, size = 3) +
    scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "#D7301F"),
                       name = NULL) +
    labs(title = paste("GSE80655", region, "MDD vs Control"),
         subtitle = sprintf("BH FDR < 0.05: %d; dashed line = nominal P=0.05; %d non-finite omitted",
                            n_sig, n_omitted),
         x = "log2 Fold Change", y = "-log10(P-value)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  save_sci_plot(p, paste0("Volcano_", region, ".tiff"))

  # FIX: MA plot with per-region average expression (was whole-dataset)
  region_vst <- vsd_matrix[, colnames(vsd_matrix) %in% regions_list[[region]], drop = FALSE]
  full_res$AveExpr <- rowMeans(region_vst)[match(full_res$GeneSymbol, rownames(region_vst))]
  ma_res <- full_res[is.finite(full_res$AveExpr) & is.finite(full_res$logFC), ]
  p_ma <- ggplot(ma_res, aes(x = AveExpr, y = logFC, color = Significant)) +
    geom_point(alpha = 0.6, size = 1.0) +
    scale_color_manual(values = c("Not Significant" = "grey60", "Significant" = "#D7301F"),
                       name = NULL) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    labs(title = paste("MA Plot (GSE80655,", region, ")"),
         x = "Average Expression (VST)", y = "log2 Fold Change") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  save_sci_plot(p_ma, paste0("MA_", region, ".tiff"))
}

# -----------------------------------------------------------------------------
# MODULE 9: MOFA EXPORT (VST matrix + metadata + feature annotation)
# -----------------------------------------------------------------------------
mofa_dir <- file.path(output_dir, "MOFA_inputs")
dir.create(mofa_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(vsd_matrix, file.path(mofa_dir, "GSE80655_expression_cleaned.rds"))
saveRDS(sample_metadata_final, file.path(mofa_dir, "GSE80655_metadata.rds"))
feature_annot <- data.frame(feature_id = rownames(vsd_matrix),
                            gene_symbol = rownames(vsd_matrix))
saveRDS(feature_annot, file.path(mofa_dir, "GSE80655_feature_annotation.rds"))

summary_rows <- lapply(regions, function(region) {
  z <- region_full_results[[region]]
  data.frame(Region = region,
             Control = sum(sample_metadata_final$Brain_Region == region &
                             sample_metadata_final$Diagnosis == "Control"),
             MDD = sum(sample_metadata_final$Brain_Region == region &
                         sample_metadata_final$Diagnosis == "MDD"),
             Genes_tested = sum(is.finite(z$P.Value)),
             Min_P = min(z$P.Value, na.rm = TRUE),
             Min_FDR = min(z$adj.P.Val, na.rm = TRUE),
             FDR_lt_0.05 = sum(z$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE))
})
write.csv(dplyr::bind_rows(summary_rows),
          file.path(output_dir, "GSE80655_analysis_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

cat("\nMOFA inputs saved to:", mofa_dir, "\n")
cat("=== GSE80655 FIXED pipeline complete ===\n")

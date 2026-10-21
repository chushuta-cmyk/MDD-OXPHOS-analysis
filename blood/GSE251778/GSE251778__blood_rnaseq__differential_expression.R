# ==============================================================================
# GSE251778_pipeline_fixed.R
# Whole-blood RNA-seq: MDD (n=80) vs control (n=89)
# Primary question: does the independently observed blood OXPHOS signal
# replicate in this cohort?  Analyses are overall and sex-stratified.
#
# Important design limitation: all male libraries are GPL20301/HiSeq and all
# female libraries are GPL24676/NovaSeq.  Sex and platform are perfectly
# confounded.  The overall ~ Sex + Group model controls their joint nuisance
# effect, but it cannot identify a platform-independent sex effect.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
dataset_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
data_dir <- file.path(dataset_dir, "data")
output_dir <- file.path(dataset_dir, "results_fixed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(data_dir, "GSE251778_Mokhtari_rna_counts.csv.gz")
metadata_file <- file.path(data_dir, "GSE251778_metadata.csv")
stopifnot(file.exists(count_file), file.exists(metadata_file))
PADJ_CUTOFF <- 0.05

cat("\n===== MODULE 1: INPUT AND DESIGN CHECKS =====\n")
count_dt <- fread(count_file, check.names = FALSE)
if (ncol(count_dt) < 3) stop("Count file has fewer than three columns.")

# GEO file: first column is a row counter, second is gene_ID, remaining columns
# are sample IDs such as 111_rna.
gene_col <- if ("gene_ID" %in% names(count_dt)) "gene_ID" else names(count_dt)[2]
gene_id <- sub("\\..*$", "", as.character(count_dt[[gene_col]]))
sample_cols <- setdiff(names(count_dt), names(count_dt)[1:2])
count_mat <- as.matrix(count_dt[, ..sample_cols])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- gene_id
colnames(count_mat) <- sub("_rna$", "", sample_cols)
if (anyDuplicated(rownames(count_mat))) {
  count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
}

meta <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
idx <- match(colnames(count_mat), meta$SampleID)
if (anyNA(idx)) stop("Count columns not found in metadata: ",
                     paste(colnames(count_mat)[is.na(idx)], collapse = ", "))
meta <- meta[idx, , drop = FALSE]
rownames(meta) <- meta$SampleID
stopifnot(identical(colnames(count_mat), rownames(meta)))
meta$Group <- factor(meta$Group, levels = c("Control", "MDD"))
meta$Sex <- factor(meta$Sex, levels = c("Male", "Female"))
meta$Platform <- factor(meta$Platform)

cat(sprintf("  Raw counts: %d genes x %d samples\n", nrow(count_mat), ncol(count_mat)))
print(with(meta, table(Group, Sex, Platform)))
if (nrow(meta) != 169 || sum(meta$Group == "MDD") != 80 ||
    sum(meta$Group == "Control") != 89) {
  stop("Expected 169 samples (80 MDD / 89 control).")
}
if (!all(as.character(meta$Platform) == ifelse(meta$Sex == "Male", "GPL20301", "GPL24676"))) {
  stop("Unexpected platform/sex mapping; re-check metadata before fitting.")
}
write.csv(meta, file.path(output_dir, "GSE251778_metadata_aligned.csv"), row.names = FALSE)

ensembl_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(count_mat),
                                        keytype = "ENSEMBL", column = "SYMBOL",
                                        multiVals = "first")

run_contrast <- function(counts, pheno, label, design_formula) {
  cat("\n===== ANALYSIS:", label, "=====\n")
  # Match Mokhtari et al. (EBioMedicine 2025): retain RNAs with a mean raw
  # count >10 within the analysis cohort.  Their covariates (age, BMI, RIN,
  # lymphocytes and polynuclear neutrophils) are not included in GEO, so this
  # local re-analysis is a minimal-covariate screen, not an exact reproduction.
  keep <- rowMeans(counts) > 10
  cat(sprintf("  Prefilter: %d / %d genes kept (mean raw count >10)\n",
              sum(keep), length(keep)))
  dds <- DESeqDataSetFromMatrix(countData = counts[keep, , drop = FALSE],
                                colData = pheno,
                                design = design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("Group", "MDD", "Control"),
                 alpha = PADJ_CUTOFF)
  out <- as.data.frame(res)
  out$ENSEMBL <- rownames(out)
  out$GeneSymbol <- unname(ensembl_symbol[out$ENSEMBL])
  out <- out[, c("ENSEMBL", "GeneSymbol", "baseMean", "log2FoldChange",
                 "lfcSE", "stat", "pvalue", "padj")]
  out <- out[order(out$pvalue), ]
  write.csv(out, file.path(output_dir, paste0("GSE251778_DEG_", label, ".csv")),
            row.names = FALSE)
  sig <- out[!is.na(out$padj) & out$padj < PADJ_CUTOFF, ]
  write.csv(sig, file.path(output_dir,
                           paste0("GSE251778_DEG_", label, "_adjP0.05.csv")),
            row.names = FALSE)
  cat(sprintf("  Tested: %d; BH FDR <0.05: %d\n", nrow(out), nrow(sig)))

  plot_df <- out[is.finite(out$log2FoldChange) & is.finite(out$pvalue) &
                   out$pvalue > 0, ]
  plot_df$Significance <- ifelse(!is.na(plot_df$padj) & plot_df$padj < PADJ_CUTOFF,
                                "FDR < 0.05", "FDR >= 0.05")
  top_lab <- head(plot_df[order(plot_df$pvalue), ], 10)
  p <- ggplot(plot_df, aes(log2FoldChange, -log10(pvalue), color = Significance)) +
    geom_point(alpha = 0.55, size = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#2C7FB8") +
    geom_text(data = top_lab, aes(label = GeneSymbol), color = "black",
              check_overlap = TRUE, vjust = -0.6, size = 3, na.rm = TRUE) +
    scale_color_manual(values = c("FDR < 0.05" = "#D7301F",
                                  "FDR >= 0.05" = "grey72")) +
    labs(title = paste("GSE251778", label, "MDD vs control"),
         subtitle = sprintf("%d genes at BH FDR < 0.05; dashed line = nominal P=0.05",
                            nrow(sig)),
         x = "log2 fold change (MDD - control)", y = expression(-log[10](italic(P))),
         color = NULL) +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  ggsave(file.path(output_dir, paste0("GSE251778_volcano_", label, ".pdf")),
         p, width = 8, height = 6)

  target <- data.frame()
  if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    ranked <- out[!is.na(out$GeneSymbol) & nzchar(out$GeneSymbol) &
                    is.finite(out$stat), c("GeneSymbol", "stat")]
    ranked <- ranked[order(abs(ranked$stat), decreasing = TRUE), ]
    ranked <- ranked[!duplicated(ranked$GeneSymbol), ]
    gene_list <- sort(setNames(ranked$stat, ranked$GeneSymbol), decreasing = TRUE)
    gsea <- suppressMessages(clusterProfiler::gseGO(
      geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
      pvalueCutoff = 1, pAdjustMethod = "BH", seed = TRUE, eps = 0,
      verbose = FALSE))
    gsea_df <- as.data.frame(gsea)
    write.csv(gsea_df,
              file.path(output_dir, paste0("GSE251778_GSEA_BP_", label, ".csv")),
              row.names = FALSE)
    target <- gsea_df[gsea_df$ID == "GO:0006119", , drop = FALSE]
    write.csv(target,
              file.path(output_dir, paste0("GSE251778_GO0006119_", label, ".csv")),
              row.names = FALSE)
    if (nrow(target)) {
      cat(sprintf("  GO:0006119 OXPHOS: NES=%.3f, P=%.3g, BH q=%.3g\n",
                  target$NES[1], target$pvalue[1], target$p.adjust[1]))
    } else {
      cat("  GO:0006119 not returned after annotation/size filters.\n")
    }
  }
  list(out = out, sig = sig, target = target)
}

overall <- run_contrast(count_mat, meta, "overall", ~ Sex + Group)
female_idx <- meta$Sex == "Female"
female <- run_contrast(count_mat[, female_idx, drop = FALSE],
                       droplevels(meta[female_idx, , drop = FALSE]),
                       "female", ~ Group)
male_idx <- meta$Sex == "Male"
male <- run_contrast(count_mat[, male_idx, drop = FALSE],
                     droplevels(meta[male_idx, , drop = FALSE]),
                     "male", ~ Group)

get_target <- function(x, field) {
  if (nrow(x$target)) x$target[[field]][1] else NA_real_
}
summary_df <- data.frame(
  Analysis = c("overall", "female", "male"),
  MDD = c(sum(meta$Group == "MDD"), sum(meta$Group[female_idx] == "MDD"),
          sum(meta$Group[male_idx] == "MDD")),
  Control = c(sum(meta$Group == "Control"), sum(meta$Group[female_idx] == "Control"),
              sum(meta$Group[male_idx] == "Control")),
  Genes_tested = c(nrow(overall$out), nrow(female$out), nrow(male$out)),
  FDR_DEGs = c(nrow(overall$sig), nrow(female$sig), nrow(male$sig)),
  OXPHOS_NES = c(get_target(overall, "NES"), get_target(female, "NES"),
                 get_target(male, "NES")),
  OXPHOS_P = c(get_target(overall, "pvalue"), get_target(female, "pvalue"),
               get_target(male, "pvalue")),
  OXPHOS_FDR = c(get_target(overall, "p.adjust"), get_target(female, "p.adjust"),
                 get_target(male, "p.adjust"))
)
write.csv(summary_df, file.path(output_dir, "GSE251778_analysis_summary.csv"),
          row.names = FALSE)

report <- paste0(
  "GSE251778 fixed local pipeline report\n",
  "====================================\n",
  "Samples: 169 (MDD=80, control=89; female=105, male=64)\n",
  "Primary model: ~ Sex + Group; contrast MDD vs control.\n",
  "Scope: minimal-covariate public-data re-analysis only. GEO does not provide\n",
  "the age, BMI, RIN, lymphocyte and neutrophil covariates used in the paper.\n",
  "Design limitation: Sex and GEO platform/library protocol are perfectly confounded;\n",
  "the model adjusts their joint nuisance effect but cannot separate them.\n",
  "Sex-stratified models are therefore reported as sensitivity analyses.\n\n",
  paste(capture.output(print(summary_df, row.names = FALSE)), collapse = "\n"), "\n")
writeLines(report, file.path(output_dir, "GSE251778_FIXED_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("\n===== GSE251778 FIXED pipeline complete =====\n")
cat("  Outputs:", output_dir, "\n")

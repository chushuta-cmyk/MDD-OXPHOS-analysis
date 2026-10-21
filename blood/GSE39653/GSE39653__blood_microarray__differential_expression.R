# ==============================================================================
# GSE39653_pipeline_fixed.R  (FIXED copy, 2026-09)
# ------------------------------------------------------------------------------
# Dataset : GSE39653 - PBMC transcriptome, 53 samples = HC 24 / BD 8 / MDD 21
#           GPL10558 Illumina HumanHT-12 V4.0
# Original: Savitz et al. 2013, Brain Behav Immun 31:161-171 (PMID 23064081)
# ------------------------------------------------------------------------------
# Fixes vs old pipeline:
#   [P0]  Bipolar (BD) samples (n=8) must be EXCLUDED -> analysis on 45 (21/24)
#   [P1]  GEO provides no age/sex covariates for this deposited matrix; the
#         model therefore uses diagnosis only and reports this limitation.
#   [P1]  Old Step1 called eBayes(fit2,...) with undefined fit2 -> fixed to fit
#   [P1]  Immune deconvolution method name must match installed package ->
#         use immunedeconv::deconvolute(method="mcp_counter") + full-feature BH
#   [P2]  Bad file name 'Step3_MechanismExploration.R .R' -> this fixed file
#   [P2]  Output scattered -> single dir: <project>/blood/GSE39653/results_fixed/
# Thresholds (user policy): adj.P/padj < 0.05; no logFC cutoff.
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(AnnotationDbi)
  library(illuminaHumanv4.db)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
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
output_dir   <- file.path(OUTPUT_ROOT, "blood/GSE39653/results_fixed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
PADJ_CUTOFF <- 0.05
mcpcounter_genes_file <- file.path(project_root, "reference/MCPcounter/genes.txt")

# small helper (must be defined before first use)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ==============================================================================
# MODULE 1: DATA PREPARATION & QC  (keep main steps of Step1_DataPrep_QC.R)
# ==============================================================================
cat("\n===== MODULE 1: DATA PREPARATION & QC =====\n")

# --- 1.1 Locate input: series matrix OR manual RData OR cleaned rds -----------
candidates <- c(
  file.path(DATA_ROOT, "Psychology_study/GSE39653_blood_data_manual.RData"),
  file.path(DATA_ROOT, "Psychology_study/GSE39653_clean_data.rds"),
  file.path(DATA_ROOT, "GSE39653/GSE39653_series_matrix.txt")
)
input_path <- candidates[file.exists(candidates)][1]
cat("  Input:", ifelse(is.na(input_path), "NOT FOUND - edit path", input_path), "\n")

expr_mat <- pheno_df <- NULL
if (!is.na(input_path) && grepl("\\.RData$", input_path)) {
  e <- new.env(); load(input_path, envir = e)
  # The supplied RData contains two top-level objects: blood_expr and
  # blood_meta.  Selecting ls(e)[1] loaded expression only and lost metadata.
  object_names <- ls(e)
  if (all(c("blood_expr", "blood_meta") %in% object_names)) {
    expr_mat <- get("blood_expr", envir = e)
    pheno_df <- get("blood_meta", envir = e)
  } else {
    dat <- get(object_names[1], envir = e)
    if (is.list(dat) && !is.data.frame(dat)) {
      expr_mat <- dat$expr %||% dat$expression %||% dat$E
      pheno_df <- dat$pheno %||% dat$metadata %||% dat$targets
    } else expr_mat <- dat
  }
} else if (!is.na(input_path) && grepl("\\.rds$", input_path)) {
  dat <- readRDS(input_path)
  if (is.list(dat) && !is.data.frame(dat)) {
    expr_mat <- dat$expr %||% dat$expression
    pheno_df <- dat$pheno %||% dat$metadata
  } else expr_mat <- dat
} else if (!is.na(input_path)) {
  # series matrix: parse manually (avoid GEOquery slowness if possible)
  gse <- tryCatch(getGEO(filename = input_path, getGPL = TRUE), error = function(e) NULL)
  if (!is.null(gse)) {
    gse <- if (is.list(gse)) gse[[1]] else gse
    expr_mat <- exprs(gse); pheno_df <- pData(gse)
  }
}
# (helper `%||%` defined in SETTINGS above)

# --- 1.2 Build diagnosis vector from available metadata ------------------------
if (is.null(expr_mat) || is.null(pheno_df)) stop("Could not load both expression and phenotype data.")
if ("geo_accession" %in% colnames(pheno_df)) rownames(pheno_df) <- pheno_df$geo_accession
pheno_df <- pheno_df[match(colnames(expr_mat), rownames(pheno_df)), , drop = FALSE]
if (anyNA(rownames(pheno_df))) stop("Expression and phenotype samples do not align.")

if (!is.null(pheno_df) && nrow(pheno_df) > 0) {
  # Try title/characteristics parsing
  if ("title" %in% colnames(pheno_df)) {
    pheno_df$Diagnosis <- ifelse(grepl("MDD|depress", pheno_df$title, ignore.case = TRUE), "MDD",
                          ifelse(grepl("Bip|BD", pheno_df$title, ignore.case = TRUE), "BD",
                          ifelse(grepl("HC|Control|Ctrl|Healthy", pheno_df$title, ignore.case = TRUE), "Control", NA)))
  }
} else {
  # Fallback: manual mapping from GEO titles (HC/MDD/BD)
  cat("  No parsed metadata; expecting user-supplied pheno in RData.\n")
}

# --- 1.3 FIX [P0]: EXCLUDE BIPOLAR (BD) samples ---------------------------------
keep <- pheno_df$Diagnosis %in% c("MDD", "Control")
expr_mdd_ctrl <- expr_mat[, keep]
pheno_mdd_ctrl <- pheno_df[keep, , drop = FALSE]
cat("  After BD exclusion:\n"); print(table(pheno_mdd_ctrl$Diagnosis))

# --- 1.4 Numeric matrix, filter low-expressed probes (keep main step) ----------
storage.mode(expr_mdd_ctrl) <- "numeric"
expr_final <- expr_mdd_ctrl[rowMeans(expr_mdd_ctrl > 0, na.rm = TRUE) >= 0.1, ]
if (median(expr_final, na.rm = TRUE) > 10) expr_final <- log2(expr_final + 1)
cat(sprintf("  Final matrix: %d probes x %d samples\n", nrow(expr_final), ncol(expr_final)))
stopifnot(identical(colnames(expr_final), rownames(pheno_mdd_ctrl)))

# Map and collapse to one probe per gene (highest mean expression) before DE.
probe_symbol <- AnnotationDbi::mapIds(
  illuminaHumanv4.db, keys = rownames(expr_final), column = "SYMBOL",
  keytype = "PROBEID", multiVals = "first")
valid <- !is.na(probe_symbol) & probe_symbol != ""
gene_expr <- expr_final[valid, , drop = FALSE]
probe_symbol <- probe_symbol[valid]
probe_mean <- rowMeans(gene_expr, na.rm = TRUE)
ord <- order(probe_symbol, -probe_mean)
gene_expr <- gene_expr[ord, , drop = FALSE]
probe_symbol <- probe_symbol[ord]
gene_expr <- gene_expr[!duplicated(probe_symbol), , drop = FALSE]
rownames(gene_expr) <- probe_symbol[!duplicated(probe_symbol)]
cat(sprintf("  Gene-level matrix: %d genes x %d samples\n", nrow(gene_expr), ncol(gene_expr)))

saveRDS(expr_final,      file.path(output_dir, "GSE39653_expr_clean.rds"))
saveRDS(gene_expr,       file.path(output_dir, "GSE39653_gene_expression.rds"))
saveRDS(pheno_mdd_ctrl,  file.path(output_dir, "GSE39653_metadata.rds"))

# ==============================================================================
# MODULE 2: DIFFERENTIAL EXPRESSION (limma MDD vs Control, keep main steps)
# ==============================================================================
cat("\n===== MODULE 2: DIFFERENTIAL EXPRESSION (limma) =====\n")

pheno_mdd_ctrl$Diagnosis <- factor(pheno_mdd_ctrl$Diagnosis, levels = c("Control", "MDD"))
design <- model.matrix(~ Diagnosis, data = pheno_mdd_ctrl)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(gene_expr, design)
fit <- eBayes(fit)                       # FIX: was eBayes(fit2,...)
deg <- topTable(fit, coef = "DiagnosisMDD", number = Inf, adjust.method = "BH")
deg$GeneSymbol <- rownames(deg)
cat(sprintf("  DEGs adj.P < 0.05: %d\n", sum(deg$adj.P.Val < PADJ_CUTOFF)))
write.csv(deg, file.path(output_dir, "GSE39653_all_DEG_results.csv"), row.names = FALSE)
write.csv(deg[deg$adj.P.Val < PADJ_CUTOFF, ],
          file.path(output_dir, "GSE39653_sig_DEGs_adjP0.05.csv"), row.names = FALSE)

plot_deg <- deg[is.finite(deg$logFC) & is.finite(deg$P.Value) & deg$P.Value > 0, ]
plot_deg$Significance <- ifelse(plot_deg$adj.P.Val < PADJ_CUTOFF,
                                "FDR < 0.05", "FDR >= 0.05")
p_volcano <- ggplot(plot_deg, aes(logFC, -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.6, size = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#2C7FB8") +
  scale_color_manual(values = c("FDR < 0.05" = "#D7301F", "FDR >= 0.05" = "grey72")) +
  labs(title = "GSE39653 PBMC: MDD vs control",
       subtitle = sprintf("%d genes at BH FDR < 0.05; dashed line = nominal P=0.05",
                          sum(deg$adj.P.Val < PADJ_CUTOFF)),
       x = "MDD - control (log2 scale)", y = expression(-log[10](italic(P))),
       color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(output_dir, "GSE39653_volcano.tiff"), p_volcano,
       width = 8, height = 6, dpi = 300, compression = "lzw")
ggsave(file.path(output_dir, "GSE39653_volcano.pdf"), p_volcano,
       width = 8, height = 6)

# ==============================================================================
# MODULE 3: IMMUNE VALIDATION (keep main steps of Step2_ImmuneValidation.R)
# ==============================================================================
cat("\n===== MODULE 3: IMMUNE VALIDATION (FIXED statistics) =====\n")

immune_expr <- gene_expr

immune_res <- NULL
if (requireNamespace("immunedeconv", quietly = TRUE)) {
  immune_res <- tryCatch(
    immunedeconv::deconvolute(immune_expr, method = "mcp_counter"),
    error = function(e) NULL)
}
if (is.null(immune_res) && requireNamespace("MCPcounter", quietly = TRUE)) {
  cat("  immunedeconv unavailable; using MCPcounter directly.\n")
  immune_res <- tryCatch({
    if (!file.exists(mcpcounter_genes_file)) stop("Local MCPcounter genes.txt is missing")
    marker_genes <- read.table(mcpcounter_genes_file, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE, check.names = FALSE,
                               colClasses = "character")
    tmp <- MCPcounter::MCPcounter.estimate(immune_expr, featuresType = "HUGO_symbols",
                                           genes = marker_genes)
    data.frame(cell_type = rownames(tmp), tmp, check.names = FALSE)
  }, error = function(e) {
    cat("  MCPcounter could not load its marker signatures:", e$message, "\n")
    NULL
  })
}
if (is.null(immune_res)) {
  cat("  Deconvolution failed; check immunedeconv / probe IDs.\n")
} else {
  immune_long <- immune_res %>%
    pivot_longer(-cell_type, names_to = "Sample", values_to = "Abundance")
  pheno_mdd_ctrl$Sample <- rownames(pheno_mdd_ctrl)
  immune_long <- immune_long %>%
    left_join(pheno_mdd_ctrl[, c("Sample", "Diagnosis")], by = "Sample") %>%
    filter(!is.na(Diagnosis))

  # FIX: BH over ALL cell types (full-feature, not post-hoc selection)
  stats <- immune_long %>%
    group_by(cell_type) %>%
    summarise(
      MDD_mean  = mean(Abundance[Diagnosis == "MDD"], na.rm = TRUE),
      Ctrl_mean = mean(Abundance[Diagnosis == "Control"], na.rm = TRUE),
      fold_change = MDD_mean / Ctrl_mean,
      p_value   = tryCatch(wilcox.test(Abundance[Diagnosis == "MDD"],
                                       Abundance[Diagnosis == "Control"])$p.value,
                           error = function(e) NA_real_),
      .groups = "drop") %>%
    mutate(padj = p.adjust(p_value, method = "BH")) %>%
    arrange(padj)
  cat("  Immune stats (full-feature BH), top 6:\n"); print(head(stats, 6))
  cat(sprintf("  Immune cell types padj < 0.05: %d\n", sum(stats$padj < PADJ_CUTOFF)))
  write.csv(stats, file.path(output_dir, "GSE39653_immune_fullBH_stats.csv"), row.names = FALSE)
}

# ==============================================================================
# MODULE 4: RANK-BASED GO-BP GSEA
# ==============================================================================
cat("\n===== MODULE 4: RANK-BASED GO-BP GSEA =====\n")
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
  ranks <- setNames(deg$t, deg$GeneSymbol)
  ranks <- ranks[is.finite(ranks) & !duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  gsea <- tryCatch(
    gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
          pAdjustMethod = "BH", pvalueCutoff = 1, seed = 1, verbose = FALSE),
    error = function(e) {
      cat("  GSEA failed:", e$message, "\n")
      NULL
    })
  if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
    gsea_df <- as.data.frame(gsea)
    write.csv(gsea_df, file.path(output_dir, "GSE39653_GSEA_BP.csv"), row.names = FALSE)
    cat(sprintf("  GO-BP terms tested: %d; BH FDR < 0.05: %d\n",
                nrow(gsea_df), sum(gsea_df$p.adjust < PADJ_CUTOFF, na.rm = TRUE)))
  }
} else {
  cat("  clusterProfiler/org.Hs.eg.db unavailable; GSEA skipped.\n")
}

# ==============================================================================
# MODULE 5: REPORT
# ==============================================================================
report <- sprintf(
"GSE39653 FIXED pipeline report (2026-09)
===========================================
Input      : 53 samples (HC24/BD8/MDD21, Savitz 2013)
FIX [P0]   : BD excluded -> %d samples (%d MDD / %d Control)
DEGs adj.P<0.05 (MDD vs Ctrl): %d
Immune cell types padj<0.05  : %d  (full-feature BH)
Note       : small n (21 vs 24) = low power; treat as validation only
",
  sum(pheno_mdd_ctrl$Diagnosis == "MDD") + sum(pheno_mdd_ctrl$Diagnosis == "Control"),
  sum(pheno_mdd_ctrl$Diagnosis == "MDD"), sum(pheno_mdd_ctrl$Diagnosis == "Control"),
  sum(deg$adj.P.Val < PADJ_CUTOFF),
  if (exists("stats")) sum(stats$padj < PADJ_CUTOFF) else NA)
writeLines(report, file.path(output_dir, "GSE39653_FIXED_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("\n===== GSE39653 FIXED pipeline complete =====\n")
cat("  Outputs:", output_dir, "\n")

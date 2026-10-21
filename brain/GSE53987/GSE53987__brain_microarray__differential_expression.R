# =============================================================================
# GSE53987 - FIXED pipeline copy
# Brain transcriptomics (MDD vs Control), Affymetrix HG-U133 Plus 2.0 (GPL570)
# Regions: Hippocampus, Pre-frontal Cortex BA46, Associative Striatum
#
# Main modules (same structure as original GSE53987_reorganized.r):
#   M0 setup -> M1 loading -> M2 log2 check -> M3 quantile norm -> M4 filtering
#   -> M5 phenotype -> M6 region DEA (limma) -> M7 annotation -> M8 gene dedup
#   -> M9 WGCNA export -> M10 MOFA export -> M11 visualization
#
# Fixes vs original:
#   1. Filtering : removed aggressive median / variance cutoffs;
#                  mild expression filter (log2 > 3 in >=10% samples)
#   2. DEG calls : adjusted P (BH, adj.P.Val < 0.05) instead of raw P < 0.01;
#                  no logFC cutoff (MDD effects are small, keep all adj-sig genes)
#   3. Dedup     : always keep probe with highest mean expression
#                  (consistent across DEG, WGCNA and MOFA exports)
#   4. Normalization : optional skip of quantile normalization when the
#                  series matrix is already RMA-processed (GEO default)
# =============================================================================

# -----------------------------------------------------------------------------
# MODULE 0: SETUP
# -----------------------------------------------------------------------------
library(stringr)
library(limma)
library(edgeR)
library(AnnotationDbi)
library(hgu133plus2.db)

# Data file candidates (edit paths as needed)
data_candidates <- c(
  "data/GSE53987_series_matrix.txt",
  "./GSE53987_series_matrix.txt",
  "/Users/tan/Developer/projects/r_project/Study/Psychology_study/datasets/GSE53987_series_matrix.txt",
  "/Users/tan/mdd_multilayer/scripts/01_brain_transcriptomics/GSE53987/data/GSE53987_series_matrix.txt"
)
find_file <- function(candidates) {
  for (p in candidates) if (file.exists(p)) return(p)
  stop("Data file not found. Tried:\n  ", paste(candidates, collapse = "\n  "))
}
file_path <- find_file(data_candidates)

project_root <- "/Users/tan/Documents/MDD-MultiOmics/transcriptomics"
output_dir <- file.path(project_root, "brain/GSE53987/results_fixed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set to FALSE if the series matrix is already RMA-normalized (GEO default).
# FIX: avoids double normalization of microarray data.
QUANTILE_NORMALIZE <- FALSE

cat("=== GSE53987 FIXED pipeline started ===\n")

# -----------------------------------------------------------------------------
# MODULE 1: DATA LOADING
# -----------------------------------------------------------------------------
raw_data <- read.delim(
  file = file_path, skip = 66, header = TRUE, sep = "\t",
  fill = TRUE, check.names = FALSE, as.is = TRUE
)
stopifnot("ID_REF" %in% colnames(raw_data))

rownames(raw_data) <- raw_data[, "ID_REF"]
expr_matrix <- as.matrix(raw_data[, -which(colnames(raw_data) == "ID_REF")])
storage.mode(expr_matrix) <- "numeric"
cat(sprintf("Loaded expression matrix: %d probes x %d samples\n",
            nrow(expr_matrix), ncol(expr_matrix)))

# -----------------------------------------------------------------------------
# MODULE 2: LOG2 TRANSFORMATION (auto-detect)
# -----------------------------------------------------------------------------
# FIX: only transform when data is clearly linear scale (max > 50).
#      RMA-processed series matrices are usually already on log2 scale.
if (max(expr_matrix, na.rm = TRUE) > 50) {
  expr_matrix <- log2(expr_matrix + 1)
  cat("Raw linear data detected -> log2(x+1) applied.\n")
} else {
  cat("Data already on log2 scale -> skipped log2.\n")
}

# -----------------------------------------------------------------------------
# MODULE 3: QUANTILE NORMALIZATION
# -----------------------------------------------------------------------------
# FIX: skip when data is already RMA-normalized to avoid double normalization.
if (QUANTILE_NORMALIZE) {
  cat("Running quantile normalization (between arrays)...\n")
  expr_matrix <- normalizeBetweenArrays(expr_matrix, method = "quantile")
} else {
  cat("Quantile normalization skipped (RMA-preprocessed data).\n")
}

# -----------------------------------------------------------------------------
# MODULE 4: EXPRESSION FILTERING
# -----------------------------------------------------------------------------
# FIX 1: only remove probes with missing values.
expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 0, ]

# FIX 2: mild low-expression filter (log2 expression > 3 in >= 10% of samples),
#        replaces the old aggressive median / variance cutoffs.
keep_probes <- rowSums(expr_matrix > 3) >= 0.1 * ncol(expr_matrix)
expr_matrix <- expr_matrix[keep_probes, ]
cat(sprintf("After low-expression filter: %d probes x %d samples\n",
            nrow(expr_matrix), ncol(expr_matrix)))

# -----------------------------------------------------------------------------
# MODULE 5: PHENOTYPE PARSING (diagnosis, region, tetrad, covariates)
# -----------------------------------------------------------------------------
geo_lines <- readLines(file_path)

# 5a. Brain region + diagnosis (from source_name, same as original)
pheno_line <- grep("^!Sample_source_name_ch1", geo_lines, value = TRUE)
if (length(pheno_line) == 0) stop("Metadata line not found.")

labels <- gsub('"', '', strsplit(pheno_line, "\t")[[1]][-1])
keep <- grepl("major depressive disorder|control", labels, ignore.case = TRUE)
selected_labels <- labels[keep]

split_labels <- str_split_fixed(selected_labels, ", ", n = 2)
brain_area <- split_labels[, 1]
disease <- split_labels[, 2]

simple <- disease == ""
brain_area[simple] <- "hippocampus"
disease[simple] <- split_labels[simple, 1]
disease_key <- tolower(trimws(gsub("_", " ", disease)))
diagnosis <- ifelse(grepl("major depressive disorder|^mdd$", disease_key),
                    "MDD", ifelse(grepl("control", disease_key), "Control", NA_character_))
if (anyNA(diagnosis)) stop("Unrecognised diagnosis label(s): ",
                           paste(unique(disease[is.na(diagnosis)]), collapse = ", "))

brain_area <- gsub(" ", "_", brain_area)
brain_area <- gsub("\\(|\\)", "", brain_area)
brain_area[grepl("Pre-frontal_cortex", brain_area)] <- "Pre-frontal_cortex_BA46"
brain_area[grepl("Associative_striatum", brain_area)] <- "Associative_striatum"
brain_area[grepl("hippocampus", brain_area)] <- "hippocampus"

# 5b. Tetrad ID from sample titles (format: "MDD_hip_10", "control_pfc_3", ...)
# FIX: the original study (Lanz 2019) used a matched tetrad design; the trailing
#      number in the sample title identifies the tetrad, which we use as a
#      blocking factor in the DEA model (recovers paired-study power).
title_line <- grep("^!Sample_title", geo_lines, value = TRUE)
all_titles <- gsub('"', '', strsplit(title_line[1], "\t")[[1]][-1])
titles_kept <- all_titles[keep]
tetrad_id <- suppressWarnings(as.numeric(sub(".*_([0-9]+)$", "\\1", titles_kept)))

# 5c. Covariates from characteristics lines (age, gender, pmi, ph, rin)
# FIX: the original study adjusted for age, sex, PMI, pH as covariates.
char_lines <- grep("^!Sample_characteristics_ch1", geo_lines, value = TRUE)
extract_chars <- function(field) {
  hit <- NULL
  for (l in char_lines) {
    vals <- gsub('"', '', strsplit(l, "\t")[[1]][-1])
    if (any(grepl(paste0("^", field, ":"), vals, ignore.case = TRUE))) { hit <- vals; break }
  }
  if (is.null(hit)) return(NULL)
  out <- sub(paste0("^", field, ":\\s*"), "", hit, ignore.case = TRUE)
  out[keep]
}
age    <- as.numeric(extract_chars("age"))
gender <- extract_chars("gender")
pmi    <- as.numeric(extract_chars("pmi"))
ph     <- as.numeric(extract_chars("ph"))
rin    <- as.numeric(extract_chars("rin"))

# 5d. Build metadata
sample_metadata <- data.frame(
  GSM_ID = colnames(expr_matrix)[keep],
  Tetrad = tetrad_id,
  Age = age, Gender = gender, PMI = pmi, pH = ph, RIN = rin,
  Brain_Area = brain_area, Disease = disease,
  Diagnosis = diagnosis,
  stringsAsFactors = FALSE
)

expr_data <- expr_matrix[, keep]
cat("Sample table (region x diagnosis):\n")
print(table(brain_area, disease))
cat("Covariate availability (non-missing):\n")
print(colSums(!is.na(sample_metadata[, c("Age", "Gender", "PMI", "pH", "RIN")])))
cat(sprintf("Unique tetrads: %d\n", length(unique(sample_metadata$Tetrad))))

regions <- c("hippocampus", "Pre-frontal_cortex_BA46", "Associative_striatum")

# Build one gene-level matrix BEFORE differential expression.  The previous
# version fitted 54,675 probes and annotated only significant probes, which
# made the complete tables unusable for cross-cohort gene-level analyses.
# Use the same deterministic rule as the other microarray pipelines: retain
# the highest-mean probe for each gene symbol.
probe_ids_all <- rownames(expr_data)
symbols_all <- suppressWarnings(
  mapIds(hgu133plus2.db, keys = probe_ids_all, column = "SYMBOL",
         keytype = "PROBEID", multiVals = "first")
)
if (sum(!is.na(symbols_all)) < 10) {
  symbols_alt <- suppressWarnings(
    mapIds(hgu133plus2.db, keys = toupper(probe_ids_all), column = "SYMBOL",
           keytype = "PROBEID", multiVals = "first")
  )
  if (sum(!is.na(symbols_alt)) > sum(!is.na(symbols_all))) symbols_all <- symbols_alt
}
keep_annot <- !is.na(symbols_all) & nzchar(symbols_all)
probe_mean <- rowMeans(expr_data[keep_annot, , drop = FALSE], na.rm = TRUE)
ord <- order(symbols_all[keep_annot], -probe_mean)
expr_unique <- expr_data[keep_annot, , drop = FALSE][ord, , drop = FALSE]
gene_symbols <- symbols_all[keep_annot][ord]
expr_unique <- expr_unique[!duplicated(gene_symbols), , drop = FALSE]
rownames(expr_unique) <- gene_symbols[!duplicated(gene_symbols)]
cat(sprintf("Gene-level matrix before DE: %d genes x %d samples\n",
            nrow(expr_unique), ncol(expr_unique)))

# -----------------------------------------------------------------------------
# MODULE 6: REGION-STRATIFIED DIFFERENTIAL EXPRESSION (limma)
# -----------------------------------------------------------------------------
deg_results <- list()
region_full_results <- list()  # complete topTable per region (for volcano plots)
region_model_metadata <- list()

for (region in regions) {
  cat("\n--- Processing region:", region, "---\n")

  gsm_ids <- sample_metadata$GSM_ID[sample_metadata$Brain_Area == region]
  region_expr <- expr_unique[, gsm_ids, drop = FALSE]
  region_meta <- sample_metadata[sample_metadata$Brain_Area == region, ]
  region_expr <- region_expr[, region_meta$GSM_ID, drop = FALSE]

  # FIX: follow the original study (Lanz 2019) - matched tetrad design + covariates.
  #      Model: Diagnosis + Age + PMI + pH + Tetrad (tetrad = blocking).
  region_meta$Diagnosis <- factor(region_meta$Diagnosis,
                                  levels = c("Control", "MDD"))
  region_meta$Tetrad_f <- factor(region_meta$Tetrad)

  # FIX (rank deficiency): some tetrads have only ONE diagnosis in this region
  #      (e.g. hippocampus: tetrad 5/18 = Control only, tetrad 9 = MDD only).
  #      Such tetrads make Diagnosis collinear with a tetrad dummy, so the
  #      design matrix is singular -> "Coefficients not estimable: DiagnosisMDD".
  #      Standard paired-analysis solution: keep only tetrads with BOTH diagnoses.
  tet_counts <- table(region_meta$Tetrad_f, region_meta$Diagnosis)
  paired_tetrads <- rownames(tet_counts)[apply(tet_counts, 1, function(r) all(r > 0))]
  n_dropped <- nrow(region_meta) - sum(region_meta$Tetrad_f %in% paired_tetrads)
  if (n_dropped > 0) {
    cat(sprintf("  Dropping %d sample(s) from single-diagnosis tetrads -> %d paired tetrads remain\n",
                n_dropped, length(paired_tetrads)))
    region_meta <- region_meta[region_meta$Tetrad_f %in% paired_tetrads, , drop = FALSE]
    region_expr <- region_expr[, region_meta$GSM_ID, drop = FALSE]
  }
  region_meta$Tetrad_f <- droplevels(region_meta$Tetrad_f)
  region_meta$Diagnosis <- droplevels(region_meta$Diagnosis)
  region_model_metadata[[region]] <- region_meta

  # NOTE: Gender is NOT included as a covariate. The tetrad design is gender-
  #      matched (MDD and Control of the same tetrad share the same sex), so the
  #      Gender column is an exact linear combination of the tetrad dummies and
  #      would make the design matrix singular ("Coefficients not estimable").
  cov_complete <- sum(complete.cases(region_meta[, c("Age", "PMI", "pH")])) == nrow(region_meta)
  if (cov_complete && length(unique(region_meta$Tetrad_f)) > 1) {
    design <- model.matrix(~ Diagnosis + Age + PMI + pH + Tetrad_f,
                           data = region_meta)
    cat("  Model: ~ Diagnosis + Age + PMI + pH + Tetrad (matched design, gender matched within tetrad)\n")
  } else {
    design <- model.matrix(~ Diagnosis, data = region_meta)
    cat("  Covariates/tetrad incomplete -> simple model: ~ Diagnosis\n")
  }
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient design in region: ", region)

  fit <- lmFit(region_expr, design)
  fit2 <- eBayes(fit)

  # FIX: use adjusted P (BH) < 0.05 as the significance threshold.
  #      No logFC cutoff (MDD effect sizes are small).
  tt_all <- topTable(fit2, coef = "DiagnosisMDD", number = Inf, adjust.method = "fdr")
  tt_all$GeneSymbol <- rownames(tt_all)
  region_full_results[[region]] <- tt_all  # keep ALL tested genes for volcano
  out_all <- tt_all
  write.csv(out_all,
            file.path(output_dir, paste0("GSE53987_all_DE_results_", region, ".csv")),
            row.names = FALSE)
  tt <- tt_all[tt_all$adj.P.Val < 0.05, ]

  if (nrow(tt) == 0) {
    cat("  No adj-significant genes in this region.\n")
    deg_results[[region]] <- data.frame()
    next
  }

  # --- MODULE 7: FORMAT SIGNIFICANT GENE-LEVEL RESULTS ---
  tt <- tt[!is.na(rownames(tt)) & !is.na(tt$adj.P.Val) & rownames(tt) != "", , drop = FALSE]
  if (nrow(tt) == 0) {
    cat("  No valid rows left after annotation filter.\n")
    deg_results[[region]] <- data.frame()
    next
  }
  res <- tt[, c("logFC", "P.Value", "adj.P.Val")]
  res$GeneSymbol <- rownames(tt)
  res$Region <- region

  deg_results[[region]] <- res
  cat(sprintf("  adj.P < 0.05 DEGs (unique genes): %d\n", nrow(res)))
}

# -----------------------------------------------------------------------------
# MODULE 8: MERGE AND SAVE DEG RESULTS
# -----------------------------------------------------------------------------
all_degs <- do.call(rbind, deg_results)
deg_file <- file.path(output_dir, "GSE53987_MDD_vs_Control_BrainRegion_DEGs.csv")
if (nrow(all_degs) > 0) {
  write.csv(all_degs, file = deg_file, row.names = FALSE)
  cat("\nDEG results saved to:", deg_file, "\n")
} else {
  write.csv(data.frame(GeneSymbol = character(),
                       logFC = numeric(), P.Value = numeric(), adj.P.Val = numeric(),
                       Region = character()), deg_file, row.names = FALSE)
  cat("\nNo adj-significant DEGs across regions; schema-valid empty CSV saved.\n")
}

# -----------------------------------------------------------------------------
# MODULE 9: GENE-LEVEL MATRIX (shared by WGCNA and MOFA exports)
# -----------------------------------------------------------------------------
# The same expr_unique object created before DE is reused for WGCNA/MOFA.
cat("\nGene-level expression matrix reused for exports:", dim(expr_unique), "\n")

# -----------------------------------------------------------------------------
# MODULE 10: WGCNA EXPORT (per region: expression + clinical)
# -----------------------------------------------------------------------------
for (region in regions) {
  region_samples <- sample_metadata$GSM_ID[sample_metadata$Brain_Area == region]
  region_expr <- expr_unique[, region_samples, drop = FALSE]
  region_meta <- sample_metadata[sample_metadata$GSM_ID %in% region_samples, ]
  region_meta <- region_meta[match(region_samples, region_meta$GSM_ID), ]

  mdd_diag <- ifelse(region_meta$Diagnosis == "MDD", 1, 0)
  clinical_df <- data.frame(
    MDD_Diagnosis = mdd_diag,
    Brain_Area = region_meta$Brain_Area,
    Disease = region_meta$Disease,
    Group = region_meta$Diagnosis,
    row.names = region_meta$GSM_ID,
    stringsAsFactors = FALSE
  )

  write.csv(region_expr, file.path(output_dir, paste0("GSE53987_", region, "_expression.csv")),
            row.names = TRUE)
  write.csv(clinical_df, file.path(output_dir, paste0("GSE53987_", region, "_clinical.csv")),
            row.names = TRUE)
  cat("WGCNA files exported for region:", region, "\n")
}

# -----------------------------------------------------------------------------
# MODULE 11: MOFA EXPORT
# -----------------------------------------------------------------------------
mofa_dir <- file.path(output_dir, "MOFA_inputs")
dir.create(mofa_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(expr_unique, file.path(mofa_dir, "GSE53987_expression_cleaned.rds"))

mofa_metadata <- sample_metadata[match(colnames(expr_unique), sample_metadata$GSM_ID), ,
                                 drop = FALSE]
rownames(mofa_metadata) <- mofa_metadata$GSM_ID
stopifnot(identical(colnames(expr_unique), rownames(mofa_metadata)))
saveRDS(mofa_metadata, file.path(mofa_dir, "GSE53987_metadata.rds"))

cat("MOFA inputs saved to:", mofa_dir, "\n")

# -----------------------------------------------------------------------------
# MODULE 12: VISUALIZATION (volcano plots, base R, no external helper)
# FIX: plot ALL tested probes per region (grey = non-significant, red = adj.P<0.05),
#      so the volcano is a full distribution even when few genes are significant.
# -----------------------------------------------------------------------------
volcano_plot <- function(deg_df, region, file) {
  if (is.null(deg_df) || nrow(deg_df) == 0) return(invisible(NULL))
  deg_df <- deg_df[is.finite(deg_df$logFC) & is.finite(deg_df$P.Value) &
                       deg_df$P.Value > 0, , drop = FALSE]
  pdf(file, width = 7, height = 6)
  plot(deg_df$logFC, -log10(deg_df$P.Value), pch = 16,
       col = ifelse(deg_df$adj.P.Val < 0.05, "red", "grey60"),
       main = paste("GSE53987", region, "MDD vs Control"),
       xlab = "log2 Fold Change", ylab = "-log10(P)")
  abline(h = -log10(0.05), lty = 2, col = "blue")
  mtext("Colour = BH FDR; dashed line = nominal P=0.05", side = 3,
        line = 0.2, cex = 0.75)
  dev.off()
}
for (region in regions) {
  full_res <- region_full_results[[region]]
  if (is.null(full_res) || nrow(full_res) == 0) next
  volcano_plot(full_res, region, file.path(output_dir, paste0("volcano_", region, ".pdf")))
}
summary_rows <- lapply(regions, function(region) {
  z <- region_full_results[[region]]
  m <- region_model_metadata[[region]]
  data.frame(Region = region, Control = sum(m$Diagnosis == "Control"),
             MDD = sum(m$Diagnosis == "MDD"),
             Genes_tested = sum(is.finite(z$P.Value)),
             Min_P = min(z$P.Value, na.rm = TRUE),
             Min_FDR = min(z$adj.P.Val, na.rm = TRUE),
             FDR_lt_0.05 = sum(z$adj.P.Val < 0.05, na.rm = TRUE))
})
write.csv(do.call(rbind, summary_rows),
          file.path(output_dir, "GSE53987_analysis_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("\n=== GSE53987 FIXED pipeline complete ===\n")

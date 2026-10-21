#!/usr/bin/env Rscript
# N2 Rebuild: GSE251778 weighted CAMERA/FRY from raw counts
# Goal: reproduce historical blood_rnaseq.tsv numbers

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# Force system path for v3.22.0
sys_libs <- .libPaths()
dev_path <- "/Users/tan/Developer/envs/R/lib/R/library"
sys_libs_clean <- sys_libs[sys_libs != dev_path]
.libPaths(sys_libs_clean)

cat("org.Hs.eg.db version:", as.character(packageVersion("org.Hs.eg.db")), "\n")

GSE_DIR <- "/Users/tan/Documents/MDD-MultiOmics/transcriptomics/blood/GSE251778"
OUT_DIR <- "/Users/tan/Documents/Codex/2026-09-19/ple/outputs/manuscript_v6_12_optionB_2026-10-13/github_repository_draft/results/cross_cohort"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Step 1: Load counts and metadata
cat("\n=== Step 1: Load data ===\n")
count_dt <- fread(file.path(GSE_DIR, "data/GSE251778_Mokhtari_rna_counts.csv.gz"), check.names = FALSE)
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

meta <- read.csv(file.path(GSE_DIR, "data/GSE251778_metadata.csv"), stringsAsFactors = FALSE, check.names = FALSE)
idx <- match(colnames(count_mat), meta$SampleID)
meta <- meta[idx, , drop = FALSE]
rownames(meta) <- meta$SampleID
stopifnot(identical(colnames(count_mat), rownames(meta)))
meta$Group <- factor(meta$Group, levels = c("Control", "MDD"))
meta$Sex <- factor(meta$Sex, levels = c("Male", "Female"))

cat("  Counts:", nrow(count_mat), "genes x", ncol(count_mat), "samples\n")
cat("  Group:", table(meta$Group), "\n")

# Step 2: Filter low-expression (mean raw count > 10, matching original pipeline)
cat("\n=== Step 2: Filter low-expression ===\n")
keep <- rowMeans(count_mat) > 10
count_mat <- count_mat[keep, ]
cat("  After mean>10 filter:", nrow(count_mat), "genes\n")

# Step 3: Map Ensembl to Symbol
cat("\n=== Step 3: Map Ensembl to Symbol ===\n")
ensembl_symbol <- mapIds(org.Hs.eg.db, keys = rownames(count_mat),
                         keytype = "ENSEMBL", column = "SYMBOL",
                         multiVals = "first")
# Remove unmapped
mapped <- !is.na(ensembl_symbol)
count_mat <- count_mat[mapped, ]
rownames(count_mat) <- toupper(ensembl_symbol[mapped])
# Collapse duplicates
if (anyDuplicated(rownames(count_mat))) {
  count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
}
cat("  After symbol mapping:", nrow(count_mat), "genes\n")

# Step 4: Get GO:0006119 gene set
cat("\n=== Step 4: GO:0006119 gene set ===\n")
go_symbols <- unique(na.omit(select(org.Hs.eg.db, keys = "GO:0006119", keytype = "GOALL", columns = "SYMBOL")$SYMBOL))
cat("  GO:0006119 symbols:", length(go_symbols), "\n")
go_symbols_upper <- toupper(go_symbols)
idx_go <- which(rownames(count_mat) %in% go_symbols_upper)
cat("  Genes in expression matrix:", length(idx_go), "\n")
cat("  Gene names:", paste(sort(rownames(count_mat)[idx_go]), collapse = ", "), "\n")

# Step 5: Design matrix
cat("\n=== Step 5: Design matrix ===\n")
# Overall model: ~ Sex + Group (matching original analysis)
design <- model.matrix(~ Sex + Group, data = meta)
cat("  Design columns:", colnames(design), "\n")
cat("  Samples:", nrow(design), "\n")

# Step 6: voom
cat("\n=== Step 6: voom ===\n")
dge <- DGEList(counts = count_mat)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = FALSE)
cat("  voom genes:", nrow(v$E), "\n")

# Step 7: CAMERA - weighted (default)
cat("\n=== Step 7: CAMERA (weighted) ===\n")
# Fixed r = 0.01
cam_fixed <- camera(v, index = idx_go, design = design, contrast = "GroupMDD", inter.gene.cor = 0.01)
cat("  Fixed r=0.01:\n")
print(cam_fixed)

# Estimated correlation
cam_est <- camera(v, index = idx_go, design = design, contrast = "GroupMDD", inter.gene.cor = NA)
cat("  Estimated correlation:\n")
print(cam_est)

# Step 8: FRY
cat("\n=== Step 8: FRY ===\n")
fry_res <- fry(v, index = idx_go, design = design, contrast = "GroupMDD")
cat("  FRY:\n")
print(fry_res)

# Step 9: Unweighted CAMERA (for comparison)
cat("\n=== Step 9: Unweighted CAMERA (for comparison) ===\n")
cam_unweighted_fixed <- camera(v, index = idx_go, design = design, contrast = "GroupMDD", inter.gene.cor = 0.01, weighted = FALSE)
cat("  Unweighted fixed r=0.01:\n")
print(cam_unweighted_fixed)

cam_unweighted_est <- camera(v, index = idx_go, design = design, contrast = "GroupMDD", inter.gene.cor = NA, weighted = FALSE)
cat("  Unweighted estimated:\n")
print(cam_unweighted_est)

fry_unweighted <- fry(v, index = idx_go, design = design, contrast = "GroupMDD")
cat("  FRY (same as above, fry is always unweighted):\n")
print(fry_unweighted)

# Step 10: Save results
cat("\n=== Step 10: Save results ===\n")
results <- data.frame(
  Test = c(
    "WEIGHTED_CAMERA_fixed_r0.01",
    "WEIGHTED_CAMERA_estimated",
    "WEIGHTED_FRY",
    "UNWEIGHTED_CAMERA_fixed_r0.01",
    "UNWEIGHTED_CAMERA_estimated",
    "UNWEIGHTED_FRY"
  ),
  NGenes = c(
    cam_fixed$NGenes,
    cam_est$NGenes,
    fry_res$NGenes,
    cam_unweighted_fixed$NGenes,
    cam_unweighted_est$NGenes,
    fry_unweighted$NGenes
  ),
  Correlation = c(
    if ("Correlation" %in% names(cam_fixed)) cam_fixed$Correlation else 0.01,
    cam_est$Correlation,
    NA,
    if ("Correlation" %in% names(cam_unweighted_fixed)) cam_unweighted_fixed$Correlation else 0.01,
    cam_unweighted_est$Correlation,
    NA
  ),
  Direction = c(
    cam_fixed$Direction,
    cam_est$Direction,
    fry_res$Direction,
    cam_unweighted_fixed$Direction,
    cam_unweighted_est$Direction,
    fry_unweighted$Direction
  ),
  PValue = c(
    cam_fixed$PValue,
    cam_est$PValue,
    fry_res$PValue,
    cam_unweighted_fixed$PValue,
    cam_unweighted_est$PValue,
    fry_unweighted$PValue
  )
)

print(results)
fwrite(results, file.path(OUT_DIR, "N2_GSE251778_weighted_camera_rebuild.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(OUT_DIR, "N2_GSE251778_sessionInfo.txt"))

cat("\n=== DONE ===\n")
cat("Output saved to:", OUT_DIR, "\n")

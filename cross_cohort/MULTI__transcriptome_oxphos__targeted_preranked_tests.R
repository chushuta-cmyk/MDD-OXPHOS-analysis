# Pre-specified cross-transcriptome check of GO:0006119 (oxidative phosphorylation).
# This is a targeted validation, not a post-hoc search across GO terms.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

repo_root <- WORK_ROOT
out_dir <- file.path(OUTPUT_ROOT, "cross_cohort")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  GSE53987_Hippocampus = "transcriptomics/brain/GSE53987/results_fixed/GSE53987_all_DE_results_hippocampus.csv",
  GSE53987_BA46 = "transcriptomics/brain/GSE53987/results_fixed/GSE53987_all_DE_results_Pre-frontal_cortex_BA46.csv",
  GSE53987_Striatum = "transcriptomics/brain/GSE53987/results_fixed/GSE53987_all_DE_results_Associative_striatum.csv",
  GSE54563_ACC = "transcriptomics/brain/GSE54563/results_fixed/GSE54563_All_DEG_Results.csv",
  GSE80655_AnCg = "transcriptomics/brain/GSE80655/results_fixed/GSE80655_all_DE_results_AnCg.csv",
  GSE80655_DLPFC = "transcriptomics/brain/GSE80655/results_fixed/GSE80655_all_DE_results_DLPFC.csv",
  GSE80655_nAcc = "transcriptomics/brain/GSE80655/results_fixed/GSE80655_all_DE_results_nAcc.csv",
  GSE98793_blood = "transcriptomics/blood/GSE98793/results_fixed/GSE98793_all_DEG_results.csv",
  GSE39653_PBMC = "transcriptomics/blood/GSE39653/results_fixed/GSE39653_all_DEG_results.csv",
  GSE19738_blood = "analysis_framework/GSE19738/results/deg_GSE19738.csv",
  GSE19738_female = "analysis_framework/GSE19738/results/deg_GSE19738_female.csv",
  GSE19738_male = "analysis_framework/GSE19738/results/deg_GSE19738_male.csv",
  GSE251778_overall = "transcriptomics/blood/GSE251778/results_fixed/GSE251778_DEG_overall.csv",
  GSE251778_female = "transcriptomics/blood/GSE251778/results_fixed/GSE251778_DEG_female.csv",
  GSE251778_male = "transcriptomics/blood/GSE251778/results_fixed/GSE251778_DEG_male.csv")
file_labels <- names(files)
files <- setNames(file.path(repo_root, unname(files)), file_labels)
if (any(!file.exists(files))) stop("Missing fixed transcriptome result(s).")

go_map <- AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0006119",
                                keytype = "GOALL", columns = "SYMBOL")
oxphos <- unique(na.omit(go_map$SYMBOL))

one <- function(file, label) {
  x <- fread(file)
  gene_col <- if ("GeneSymbol" %in% names(x)) "GeneSymbol" else
    if ("gene" %in% names(x)) "gene" else NA_character_
  if (is.na(gene_col)) stop(label, ": gene-symbol column missing")
  if ("stat" %in% names(x)) {
    score <- x$stat
  } else if ("t" %in% names(x)) {
    score <- x$t
  } else {
    score <- sign(x$logFC) * qnorm(pmax(x$P.Value, 1e-300) / 2,
                                   lower.tail = FALSE)
  }
  ranks <- setNames(score, x[[gene_col]])
  ranks <- ranks[is.finite(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))]
  ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)
  set.seed(20260919)
  y <- fgseaMultilevel(pathways = list(Oxidative_phosphorylation_GO_0006119 = oxphos),
                       stats = ranks, minSize = 10, maxSize = 500,
                       eps = 0)
  data.table(Dataset = label, Tested_genes = length(ranks),
             SetSize = y$size, NES = y$NES, P = y$pval,
             LeadingEdge = paste(y$leadingEdge[[1]], collapse = ";"))
}

res <- rbindlist(Map(one, files, names(files)))
res[, FDR_across_datasets := p.adjust(P, method = "BH")]
fwrite(res, file.path(out_dir, "predefined_OXPHOS_across_transcriptomes.tsv"), sep = "\t")
print(res[, .(Dataset, Tested_genes, SetSize, NES, P, FDR_across_datasets)])

# Cross-validation evidence: quantify overlap among the three independent blood
# cohorts showing negative enrichment.  GSE251778 is kept in the full table but
# excluded from this core because its public metadata lack the five covariates
# used in the paper and its male/female directions differ.
le_sets <- setNames(strsplit(res$LeadingEdge, ";", fixed = TRUE), res$Dataset)
core_ids <- c("GSE98793_blood", "GSE39653_PBMC", "GSE19738_blood")
core_genes <- sort(Reduce(intersect, le_sets[core_ids]))
fwrite(data.table(GeneSymbol = core_genes),
       file.path(out_dir, "blood_OXPHOS_three_cohort_core_genes.tsv"), sep = "\t")

blood_ids <- c(core_ids, "GSE19738_female", "GSE19738_male",
               "GSE251778_female", "GSE251778_male")
pairs <- combn(blood_ids, 2, simplify = FALSE)
pairwise <- rbindlist(lapply(pairs, function(z) {
  a <- le_sets[[z[1]]]; b <- le_sets[[z[2]]]
  data.table(Dataset_1 = z[1], Dataset_2 = z[2],
             N_1 = length(a), N_2 = length(b), Overlap = length(intersect(a, b)),
             Jaccard = length(intersect(a, b)) / length(union(a, b)),
             Same_NES_direction = sign(res$NES[match(z[1], res$Dataset)]) ==
               sign(res$NES[match(z[2], res$Dataset)]))
}))
fwrite(pairwise, file.path(out_dir, "blood_OXPHOS_leading_edge_pairwise.tsv"), sep = "\t")
cat(sprintf("\nThree-cohort negative blood OXPHOS core: %d genes\n", length(core_genes)))
if (length(core_genes)) cat("  ", paste(core_genes, collapse = ", "), "\n")

core_detail <- rbindlist(lapply(blood_ids, function(id) {
  x <- fread(files[[id]])
  gene_col <- if ("GeneSymbol" %in% names(x)) "GeneSymbol" else "gene"
  stat_value <- if ("stat" %in% names(x)) x$stat else if ("t" %in% names(x)) x$t else
    sign(x$logFC) * qnorm(pmax(x$P.Value, 1e-300) / 2, lower.tail = FALSE)
  p_value <- if ("P.Value" %in% names(x)) x$P.Value else x$pvalue
  logfc <- if ("logFC" %in% names(x)) x$logFC else x$log2FoldChange
  z <- data.table(Dataset = id, GeneSymbol = x[[gene_col]],
                  logFC = logfc, Statistic = stat_value, P = p_value)
  z[GeneSymbol %in% core_genes][, Direction := ifelse(logFC > 0, "UP", "DOWN")]
}))
setorder(core_detail, GeneSymbol, Dataset)
fwrite(core_detail, file.path(out_dir, "blood_OXPHOS_core_gene_cross_dataset.tsv"),
       sep = "\t")

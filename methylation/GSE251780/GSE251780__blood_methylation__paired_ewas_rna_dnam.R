#!/usr/bin/env Rscript

# GSE251780 (EPIC methylation) paired with GSE251778 (RNA-seq).
# Processed-data analysis: MDD DMPs, pre-specified GO:0006119 promoter
# methylation, and within-person methylation--RNA concordance.
# Run from repository root:
#   Rscript analysis_framework/15_GSE251780_paired_methylation_RNA.R

# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

repo_root <- WORK_ROOT
source(file.path(repo_root, "analysis_framework", "_libpaths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(org.Hs.eg.db)
  library(missMethyl)
  library(ggplot2)
  library(ggrepel)
  library(sandwich)
  library(lmtest)
})

set.seed(251780)
base_dir <- file.path(repo_root, "analysis_framework", "GSE251780")
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(OUTPUT_ROOT, "methylation/GSE251780/results")
fig_dir <- file.path(OUTPUT_ROOT, "methylation/GSE251780/figures")
log_dir <- file.path(OUTPUT_ROOT, "methylation/GSE251780/logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(data_dir,
  "GSE251780_Mokhtari_MatrixProcessed_AverageAndPvalue.txt.gz")
metadata_file <- file.path(data_dir, "GSE251780_metadata.csv")
rna_score_file <- file.path(repo_root, "analysis_framework", "results_strict",
                             "blood_OXPHOS_sample_scores.tsv")
rna_count_file <- file.path(repo_root, "transcriptomics", "blood", "GSE251778",
                             "data", "GSE251778_Mokhtari_rna_counts.csv.gz")
stopifnot(file.exists(matrix_file), file.exists(metadata_file),
          file.exists(rna_score_file), file.exists(rna_count_file))

cat("===== GSE251780 paired methylation--RNA analysis =====\n")
gzcmd <- paste("gzip -cd", shQuote(matrix_file))
hdr <- names(fread(cmd = gzcmd, nrows = 0, check.names = FALSE))
id_col <- hdr[1]
det_cols <- grep("\\.Detection\\.Pval$", hdr, value = TRUE)
beta_cols <- setdiff(hdr[-1], det_cols)
det_ids <- sub("\\.Detection\\.Pval$", "", det_cols)
if (!setequal(det_ids, beta_cols)) stop("Detection-P and beta sample names differ")
cat("Matrix header:", length(beta_cols), "beta samples and",
    length(det_cols), "detection-P samples\n")

# Read detection P values first so failed probes can be discarded before the
# much larger beta-value object is retained in memory.
det_dt <- fread(cmd = gzcmd, select = c(id_col, det_cols),
                showProgress = interactive(), check.names = FALSE)
probe_id <- as.character(det_dt[[id_col]])
det <- as.matrix(det_dt[, ..det_cols])
storage.mode(det) <- "double"
colnames(det) <- det_ids
probe_det_pass <- rowMeans(det < 0.01, na.rm = TRUE) >= 0.95
sample_qc <- data.table(
  SampleID = det_ids,
  Failed_probe_fraction_P001 = colMeans(det >= 0.01 | !is.finite(det), na.rm = TRUE))
fwrite(sample_qc, file.path(out_dir, "GSE251780_detectionP_sample_QC.tsv"), sep = "\t")
rm(det_dt, det); gc()

ann <- as.data.frame(minfi::getAnnotation(
  IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
ann$Probe <- rownames(ann)
ann_idx <- match(probe_id, ann$Probe)
autosomal <- !is.na(ann_idx) & ann$chr[ann_idx] %in% paste0("chr", 1:22)
common_snp <- (!is.na(ann$CpG_maf[ann_idx]) & ann$CpG_maf[ann_idx] >= 0.05) |
              (!is.na(ann$SBE_maf[ann_idx]) & ann$SBE_maf[ann_idx] >= 0.05)
common_snp[is.na(common_snp)] <- FALSE
keep_probe <- probe_det_pass & autosomal & !common_snp
keep_ids <- probe_id[keep_probe]
cat("Probe QC:", sum(probe_det_pass), "detection-P pass;",
    length(keep_ids), "autosomal/non-common-SNP retained\n")
rm(probe_det_pass, autosomal, common_snp, ann_idx); gc()

beta_dt <- fread(cmd = gzcmd, select = c(id_col, beta_cols),
                 showProgress = interactive(), check.names = FALSE)
beta_probe <- as.character(beta_dt[[id_col]])
ii <- match(keep_ids, beta_probe)
if (anyNA(ii)) stop("Retained probe IDs missing from beta matrix")
beta <- as.matrix(beta_dt[ii, ..beta_cols])
storage.mode(beta) <- "double"
rownames(beta) <- keep_ids
colnames(beta) <- beta_cols
rm(beta_dt, beta_probe, ii, probe_id, keep_ids); gc()
beta[beta <= 0] <- 1e-6
beta[beta >= 1] <- 1 - 1e-6
M <- log2(beta / (1 - beta))
rm(beta); gc()

meta <- fread(metadata_file)
meta[, SampleID := as.character(SampleID)]
missing_meta <- setdiff(colnames(M), meta$SampleID)
if (length(missing_meta)) stop("Methylation samples absent from metadata: ",
                               paste(missing_meta, collapse = ", "))
meta <- meta[match(colnames(M), SampleID)]
meta[, Group := factor(Group, levels = c("Control", "MDD"))]
meta[, Sex := factor(Sex, levels = c("Female", "Male"))]
stopifnot(!anyNA(meta$Group), !anyNA(meta$Sex))

design <- model.matrix(~ Group + Sex, data = meta)
fit <- eBayes(lmFit(M, design), robust = TRUE, trend = TRUE)
dmp <- topTable(fit, coef = "GroupMDD", number = Inf, sort.by = "P")
dmp$Probe <- rownames(dmp)
ann_sub <- ann[match(dmp$Probe, ann$Probe), , drop = FALSE]
for (nm in c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group",
             "Relation_to_Island")) dmp[[nm]] <- ann_sub[[nm]]
dmp <- dmp[order(dmp$P.Value), ]
fwrite(dmp, file.path(out_dir, "GSE251780_MDD_DMP_all.tsv"), sep = "\t")
fwrite(dmp[dmp$adj.P.Val < 0.05, ],
       file.path(out_dir, "GSE251780_MDD_DMP_FDR05.tsv"), sep = "\t")

# Pre-defined oxidative phosphorylation set.
go_tab <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = "GO:0006119", keytype = "GOALL",
  columns = c("SYMBOL", "ENTREZID", "ONTOLOGYALL")))
go_tab <- go_tab[go_tab$ONTOLOGYALL == "BP", ]
ox_genes <- sort(unique(na.omit(go_tab$SYMBOL)))
ox_entrez <- unique(na.omit(as.character(go_tab$ENTREZID)))
ann_m <- ann[match(rownames(M), ann$Probe), , drop = FALSE]
gene_strings <- as.character(ann_m$UCSC_RefGene_Name)
gene_split <- strsplit(gene_strings, ";", fixed = TRUE)
is_ox <- vapply(gene_split, function(z) any(z %in% ox_genes), logical(1))
is_promoter <- grepl("TSS1500|TSS200|1stExon",
                     as.character(ann_m$UCSC_RefGene_Group))
ox_probe_idx <- which(is_ox & is_promoter)
ox_probes <- rownames(M)[ox_probe_idx]
writeLines(ox_probes, file.path(out_dir, "GSE251780_OXPHOS_promoter_probes.txt"))

# Probe-number-bias-corrected competitive methylation gene-set tests.
gs_res <- rbindlist(lapply(c(0.001, 0.005, 0.01), function(cutoff) {
  sig <- dmp$Probe[dmp$P.Value < cutoff]
  if (length(sig) < 2) return(data.table(
    CpG_P_cutoff = cutoff, N_selected_CpGs = length(sig),
    Status = "Too few CpGs"))
  z <- tryCatch(missMethyl::gsameth(
    sig.cpg = sig, all.cpg = rownames(M),
    collection = list(GO_0006119_OXPHOS = ox_entrez),
    array.type = "EPIC", plot.bias = FALSE), error = function(e) e)
  if (inherits(z, "error")) return(data.table(
    CpG_P_cutoff = cutoff, N_selected_CpGs = length(sig),
    Status = conditionMessage(z)))
  ans <- as.data.table(z, keep.rownames = "GeneSet")
  ans[, `:=`(CpG_P_cutoff = cutoff, N_selected_CpGs = length(sig), Status = "OK")]
  ans
}), fill = TRUE)
fwrite(gs_res, file.path(out_dir, "GSE251780_MDD_OXPHOS_gsameth.tsv"), sep = "\t")

# Sample-level OXPHOS promoter methylation score.  This is a targeted,
# interpretable sensitivity score; promoter methylation is not assumed to have
# a uniform functional direction for every probe.
ox_M <- M[ox_probe_idx, , drop = FALSE]
ox_z <- t(scale(t(ox_M)))
dnam_score <- as.numeric(scale(colMeans(ox_z, na.rm = TRUE)))
score_dat <- data.frame(meta, DNAm_OXPHOS_promoter_z = dnam_score)
score_fit <- lm(DNAm_OXPHOS_promoter_z ~ Group + Sex, data = score_dat)
score_ct <- lmtest::coeftest(score_fit, vcov. = sandwich::vcovHC(score_fit, type = "HC3"))
b <- unname(score_ct["GroupMDD", 1]); se <- unname(score_ct["GroupMDD", 2])
score_result <- data.table(
  Analysis = "MDD effect on OXPHOS promoter methylation mean-z score",
  N = nobs(score_fit), N_CpGs = nrow(ox_M), Beta_MDD = b, SE_HC3 = se,
  CI_low = b - 1.96 * se, CI_high = b + 1.96 * se,
  P = unname(score_ct["GroupMDD", 4]),
  Caveat = "Processed beta values; no public cell-fraction covariates")
fwrite(score_result, file.path(out_dir, "GSE251780_OXPHOS_DNAm_score_model.tsv"), sep = "\t")

# Join the already validated GSE251778 transcriptomic sample scores.
rna_scores <- fread(rna_score_file)[Cohort == "GSE251778"]
rna_scores[, SampleID := as.character(Sample)]
paired <- merge(as.data.table(score_dat),
                rna_scores[, .(SampleID, RNA_OXPHOS_all = Score_all,
                                RNA_OXPHOS_nuclear = Score_nuclear,
                                RNA_OXPHOS_core13 = Score_core13)],
                by = "SampleID", all = FALSE)
if (nrow(paired) != nrow(meta)) warning("Not every methylation sample has an RNA score")
paired_fit <- lm(RNA_OXPHOS_all ~ DNAm_OXPHOS_promoter_z + Group + Sex,
                 data = paired)
paired_ct <- lmtest::coeftest(paired_fit,
  vcov. = sandwich::vcovHC(paired_fit, type = "HC3"))
pb <- unname(paired_ct["DNAm_OXPHOS_promoter_z", 1])
pse <- unname(paired_ct["DNAm_OXPHOS_promoter_z", 2])
paired_result <- data.table(
  Analysis = "Paired RNA OXPHOS score vs promoter DNAm score",
  N = nobs(paired_fit), Beta_DNAm = pb, SE_HC3 = pse,
  CI_low = pb - 1.96 * pse, CI_high = pb + 1.96 * pse,
  P = unname(paired_ct["DNAm_OXPHOS_promoter_z", 4]),
  Expected_direction = "Negative if aggregate promoter methylation suppresses RNA")
fwrite(paired_result, file.path(out_dir, "GSE251780_paired_pathway_association.tsv"), sep = "\t")
fwrite(paired, file.path(out_dir, "GSE251780_paired_sample_scores.tsv"), sep = "\t")

# Gene-resolved paired analysis for OXPHOS genes.
counts_dt <- fread(cmd = paste("gzip -cd", shQuote(rna_count_file)),
                   check.names = FALSE)
ens <- sub("\\..*$", "", as.character(counts_dt$gene_ID))
count_cols <- grep("_rna$", names(counts_dt), value = TRUE)
counts <- as.matrix(counts_dt[, ..count_cols]); storage.mode(counts) <- "double"
rownames(counts) <- ens
colnames(counts) <- sub("_rna$", "", count_cols)
keep_rna <- rowSums(counts >= 10, na.rm = TRUE) >= 10
lcpm <- edgeR::cpm(counts[keep_rna, , drop = FALSE], log = TRUE, prior.count = 2)
ens_keep <- rownames(lcpm)
ens_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = unique(ens_keep), keytype = "ENSEMBL", columns = "SYMBOL"))
ens_map <- ens_map[!is.na(ens_map$SYMBOL) & ens_map$SYMBOL %in% ox_genes, ]
ens_map <- ens_map[match(ens_keep, ens_map$ENSEMBL, nomatch = 0L), ]
rna_symbols <- unique(ens_map$SYMBOL)
rna_gene <- do.call(rbind, lapply(rna_symbols, function(g) {
  ids <- unique(ens_map$ENSEMBL[ens_map$SYMBOL == g])
  colMeans(lcpm[intersect(ids, rownames(lcpm)), , drop = FALSE], na.rm = TRUE)
}))
rownames(rna_gene) <- rna_symbols

# A promoter probe may map to multiple genes; each listed OXPHOS gene receives
# the probe.  Duplicate probe assignments within a gene are removed.
prom_map <- rbindlist(lapply(ox_probe_idx, function(i) {
  gs <- intersect(unique(gene_split[[i]]), ox_genes)
  data.table(Probe = rownames(M)[i], Gene = gs)
}))
prom_map <- unique(prom_map)
dnam_genes <- unique(prom_map$Gene)
dnam_gene <- do.call(rbind, lapply(dnam_genes, function(g) {
  pp <- prom_map$Probe[prom_map$Gene == g]
  colMeans(M[pp, , drop = FALSE], na.rm = TRUE)
}))
rownames(dnam_gene) <- dnam_genes

common_gene <- intersect(rownames(rna_gene), rownames(dnam_gene))
common_sample <- Reduce(intersect, list(colnames(rna_gene), colnames(dnam_gene),
                                        paired$SampleID))
rna_g <- rna_gene[common_gene, common_sample, drop = FALSE]
dnam_g <- dnam_gene[common_gene, common_sample, drop = FALSE]
covdat <- as.data.frame(meta[match(common_sample, SampleID), .(Group, Sex)])

gene_res <- rbindlist(lapply(seq_along(common_gene), function(i) {
  dd <- data.frame(RNA = as.numeric(rna_g[i, ]),
                   DNAm = as.numeric(dnam_g[i, ]), covdat)
  ff <- lm(RNA ~ DNAm + Group + Sex, data = dd)
  cc <- summary(ff)$coefficients["DNAm", ]
  data.table(Gene = common_gene[i], N = nobs(ff),
             N_promoter_CpGs = sum(prom_map$Gene == common_gene[i]),
             Beta_DNAm = cc[1], SE = cc[2], T = cc[3], P = cc[4])
}))
gene_res[, FDR := p.adjust(P, "BH")]
setorder(gene_res, P)
fwrite(gene_res, file.path(out_dir, "GSE251780_paired_OXPHOS_gene_associations.tsv"),
       sep = "\t")

# Gene-set-level paired concordance with sample-label permutations within the
# four diagnosis-by-sex strata. This preserves major phenotype structure and
# accounts for correlation among OXPHOS genes.
cov_design <- model.matrix(~ Group + Sex, covdat)
rna_covfit <- lmFit(rna_g, cov_design)
dnam_covfit <- lmFit(dnam_g, cov_design)
rna_res <- rna_g - rna_covfit$coefficients %*% t(cov_design)
dnam_res <- dnam_g - dnam_covfit$coefficients %*% t(cov_design)
row_standardize <- function(x) {
  x <- x - rowMeans(x, na.rm = TRUE)
  x / sqrt(rowSums(x^2, na.rm = TRUE))
}
rz <- row_standardize(rna_res); dz <- row_standardize(dnam_res)
obs_cor <- rowSums(rz * dz, na.rm = TRUE)
obs_mean_z <- mean(atanh(pmin(pmax(obs_cor, -0.999999), 0.999999)))
strata <- split(seq_len(ncol(rz)), interaction(covdat$Group, covdat$Sex, drop = TRUE))
nperm <- 2000L
perm_stat <- numeric(nperm)
for (k in seq_len(nperm)) {
  perm <- seq_len(ncol(rz))
  for (s in strata) perm[s] <- sample(s, length(s), replace = FALSE)
  rr <- rowSums(rz[, perm, drop = FALSE] * dz, na.rm = TRUE)
  perm_stat[k] <- mean(atanh(pmin(pmax(rr, -0.999999), 0.999999)))
}
perm_result <- data.table(
  Analysis = "Mean partial within-gene RNA-promoter-DNAm correlation",
  N_pairs = length(common_sample), N_genes = length(common_gene),
  Mean_Fisher_z = obs_mean_z,
  Mean_partial_r = tanh(obs_mean_z),
  P_two_sided_permutation = (1 + sum(abs(perm_stat) >= abs(obs_mean_z))) / (nperm + 1),
  P_expected_negative = (1 + sum(perm_stat <= obs_mean_z)) / (nperm + 1),
  N_permutations = nperm)
fwrite(perm_result, file.path(out_dir, "GSE251780_paired_gene_set_permutation.tsv"),
       sep = "\t")

# Figures.
p_group <- ggplot(score_dat, aes(Group, DNAm_OXPHOS_promoter_z, fill = Group)) +
  geom_violin(alpha = 0.35, trim = FALSE) +
  geom_boxplot(width = 0.20, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.08, alpha = 0.55, size = 1.4) +
  scale_fill_manual(values = c(Control = "#3B82B7", MDD = "#C74B50")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(title = "GSE251780: OXPHOS promoter methylation in MDD",
       subtitle = "Inference adjusted for sex; processed-data sensitivity",
       x = NULL, y = "OXPHOS promoter methylation score (z)")
ggsave(file.path(fig_dir, "GSE251780_OXPHOS_DNAm_by_group.png"), p_group,
       width = 7.2, height = 6, dpi = 300)

p_pair <- ggplot(paired, aes(DNAm_OXPHOS_promoter_z, RNA_OXPHOS_all,
                             color = Group, shape = Sex)) +
  geom_point(alpha = 0.72, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black") +
  scale_color_manual(values = c(Control = "#3B82B7", MDD = "#C74B50")) +
  theme_bw(base_size = 12) +
  labs(title = "Paired blood methylation and RNA OXPHOS scores",
       subtitle = "Line is unadjusted; reported inference adjusts for diagnosis and sex",
       x = "Promoter DNAm OXPHOS score (z)", y = "RNA OXPHOS score (z)")
ggsave(file.path(fig_dir, "GSE251780_paired_DNAm_RNA_scatter.png"), p_pair,
       width = 8, height = 6.2, dpi = 300)

vol <- as.data.table(dmp)
vol[, minus_log10_p := -log10(pmax(P.Value, .Machine$double.xmin))]
vol[, Label := ifelse(seq_len(.N) <= 10, Probe, NA_character_)]
p_vol <- ggplot(vol, aes(logFC, minus_log10_p)) +
  geom_point(alpha = 0.28, size = 0.55, color = "grey55") +
  geom_point(data = vol[adj.P.Val < 0.05], color = "#C0392B", size = 0.8) +
  geom_text_repel(aes(label = Label), na.rm = TRUE, max.overlaps = Inf, size = 3) +
  theme_bw(base_size = 12) +
  labs(title = "GSE251780 MDD methylome", x = "MDD - control (M-value)",
       y = expression(-log[10](P)))
ggsave(file.path(fig_dir, "GSE251780_MDD_volcano.png"), p_vol,
       width = 8.2, height = 6.3, dpi = 300)

p_cor <- ggplot(gene_res, aes(Beta_DNAm)) +
  geom_histogram(bins = 28, fill = "#6A51A3", color = "white") +
  geom_vline(xintercept = 0, linetype = 2) + theme_bw(base_size = 12) +
  labs(title = "Paired OXPHOS gene-level methylation--RNA effects",
       subtitle = "RNA logCPM ~ promoter methylation + diagnosis + sex",
       x = "Adjusted RNA effect per promoter M-value", y = "Genes")
ggsave(file.path(fig_dir, "GSE251780_paired_gene_effect_distribution.png"), p_cor,
       width = 7.5, height = 5.8, dpi = 300)

summary_lines <- c(
  "GSE251780/GSE251778 paired processed-data analysis",
  "=================================================",
  sprintf("Methylation samples: %d (controls %d; MDD %d)", nrow(meta),
          sum(meta$Group == "Control"), sum(meta$Group == "MDD")),
  sprintf("QC-retained autosomal probes: %d", nrow(M)),
  sprintf("FDR-significant MDD CpGs: %d", sum(dmp$adj.P.Val < 0.05, na.rm = TRUE)),
  sprintf("OXPHOS promoter score CpGs: %d", nrow(ox_M)),
  sprintf("Paired methylation--RNA samples: %d", nrow(paired)),
  sprintf("Gene-resolved paired OXPHOS genes: %d", length(common_gene)),
  "Covariates: sex (age and blood-cell fractions are unavailable publicly)",
  "Interpretation must distinguish within-cohort multi-omic concordance from external replication.")
writeLines(summary_lines, file.path(out_dir, "GSE251780_analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "GSE251780_sessionInfo.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

#!/usr/bin/env Rscript

# GSE298549: processed-data sensitivity analysis of weighted life events (wLE)
# in the MDD arm of GSE201016.  This is deliberately labelled a sensitivity
# analysis because the public processed matrix does not contain estimated blood
# cell fractions.  A raw-IDAT server workflow is provided separately.
#
# Run from the repository root:
#   Rscript scripts/methylation/GSE298549__blood_methylation__weighted_life_events.R
# Optional DMR step (slower):
#   RUN_DMR=true Rscript scripts/methylation/GSE298549__blood_methylation__weighted_life_events.R

repo_root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo_root, "analysis_framework", "_libpaths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(org.Hs.eg.db)
  library(missMethyl)
  library(ggplot2)
  library(ggrepel)
  library(sandwich)
  library(lmtest)
})

set.seed(298549)
base_dir <- file.path(repo_root, "analysis_framework", "GSE201016")
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "figures")
log_dir <- file.path(base_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(data_dir, "GSE298549_mval_MDD_cases.txt.gz")
metadata_file <- file.path(data_dir, "GSE298549_metadata.csv")
if (!file.exists(matrix_file)) {
  stop("Missing processed M-value matrix: ", matrix_file,
       "\nRun analysis_framework/download_validation_datasets.sh first.")
}
if (!file.exists(metadata_file)) stop("Missing metadata: ", metadata_file)

cat("===== GSE298549 processed M-value sensitivity analysis =====\n")
cat("Reading matrix (this can require several GB of RAM)\n")
tab <- fread(cmd = paste("gzip -cd", shQuote(matrix_file)),
             showProgress = interactive())
probe_col <- names(tab)[1]
probe_id <- as.character(tab[[probe_col]])
tab[, (probe_col) := NULL]
M <- as.matrix(tab)
storage.mode(M) <- "double"
rownames(M) <- probe_id
rm(tab)
cat("Loaded:", nrow(M), "probes x", ncol(M), "samples\n")

meta <- fread(metadata_file)
meta[, ArrayID := as.character(ArrayID)]
missing_samples <- setdiff(meta$ArrayID, colnames(M))
if (length(missing_samples)) {
  stop("Metadata samples absent from matrix: ", paste(missing_samples, collapse = ", "))
}
meta <- meta[match(colnames(M), ArrayID)]
stopifnot(identical(meta$ArrayID, colnames(M)))
meta[, Sex := factor(Sex, levels = c("Female", "Male"))]
meta[, Slide := factor(sub("_.*", "", ArrayID))]
meta[, Position := factor(sub("^[^_]+_", "", ArrayID))]
meta[, Age_z := as.numeric(scale(Age))]
meta[, WLE_z := as.numeric(scale(WLE))]

complete <- complete.cases(meta[, .(WLE_z, Age_z, Sex, Slide, Position)])
meta_cc <- droplevels(as.data.frame(meta[complete]))
M_cc <- M[, complete, drop = FALSE]
cat("Complete wLE model N:", ncol(M_cc), "of", ncol(M), "\n")
if (ncol(M_cc) < 50) stop("Too few complete samples for the planned model")

ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- as.data.frame(ann)
ann$Probe <- rownames(ann)
keep_probe <- intersect(rownames(M_cc), ann$Probe[ann$chr %in% paste0("chr", 1:22)])
M_cc <- M_cc[keep_probe, , drop = FALSE]
ann <- ann[match(rownames(M_cc), ann$Probe), , drop = FALSE]
finite_fraction <- rowMeans(is.finite(M_cc))
keep_finite <- finite_fraction >= 0.95
M_cc <- M_cc[keep_finite, , drop = FALSE]
ann <- ann[keep_finite, , drop = FALSE]

# Position is estimable without exhausting the limited residual degrees of
# freedom.  Slide has 28 sparse levels among 74 complete samples, so including
# it as a fixed effect creates a 42-column design and unstable HC3 estimates.
# We therefore use age, sex and array position for the primary processed-data
# sensitivity analysis.  The sparse-slide model is retained below as an
# explicitly over-adjusted gene-level sensitivity check.
design <- model.matrix(~ WLE_z + Age_z + Sex + Position, data = meta_cc)
nonestimable <- limma::nonEstimable(design)
if (length(nonestimable)) {
  if ("WLE_z" %in% nonestimable) stop("The wLE coefficient is non-estimable")
  design <- design[, !colnames(design) %in% nonestimable, drop = FALSE]
}
if (!"WLE_z" %in% colnames(design)) stop("WLE_z coefficient missing from design")
cat("Design:", nrow(design), "samples x", ncol(design), "coefficients\n")

fit <- eBayes(lmFit(M_cc, design), robust = TRUE, trend = TRUE)
ewas <- topTable(fit, coef = "WLE_z", number = Inf, sort.by = "P")
ewas$Probe <- rownames(ewas)
anno_cols <- intersect(c("Probe", "chr", "pos", "UCSC_RefGene_Name",
                         "UCSC_RefGene_Group", "Relation_to_Island"), colnames(ann))
ewas <- merge(ewas, ann[, anno_cols, drop = FALSE], by = "Probe", all.x = TRUE,
              sort = FALSE)
ewas <- ewas[order(ewas$P.Value), ]
fwrite(ewas, file.path(out_dir, "GSE298549_wLE_EWAS_all.tsv"), sep = "\t")
fwrite(ewas[ewas$adj.P.Val < 0.05, ],
       file.path(out_dir, "GSE298549_wLE_EWAS_FDR05.tsv"), sep = "\t")

design_slide <- model.matrix(~ WLE_z + Age_z + Sex + Slide + Position,
                             data = meta_cc)
ne_slide <- limma::nonEstimable(design_slide)
if (length(ne_slide)) {
  if ("WLE_z" %in% ne_slide) stop("WLE_z is non-estimable in slide sensitivity")
  design_slide <- design_slide[, !colnames(design_slide) %in% ne_slide, drop = FALSE]
}
fit_slide <- eBayes(lmFit(M_cc, design_slide), robust = TRUE, trend = TRUE)
ewas_slide <- topTable(fit_slide, coef = "WLE_z", number = Inf, sort.by = "P")
ewas_slide$Probe <- rownames(ewas_slide)
fwrite(ewas_slide,
       file.path(out_dir, "GSE298549_wLE_EWAS_sparse_slide_sensitivity.tsv"),
       sep = "\t")

# Direct replication check for the two CpGs highlighted in the parent study.
reported <- c("cg03341655", "cg09738429")
replication <- ewas[match(reported, ewas$Probe), ]
replication$Reported_interpretation <- c(
  "GFOD2 probe reported to replicate in the independent MDD sample",
  "Near PYROXD1; reported PD discovery signal")
fwrite(replication, file.path(out_dir, "GSE298549_reported_CpG_replication.tsv"), sep = "\t")

# Pre-defined GO:0006119 gene set and probe mapping.
go_tab <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = "GO:0006119", keytype = "GOALL",
  columns = c("SYMBOL", "ONTOLOGYALL")))
ox_genes <- sort(unique(go_tab$SYMBOL[go_tab$ONTOLOGYALL == "BP" &
                                       !is.na(go_tab$SYMBOL)]))
split_genes <- strsplit(as.character(ann$UCSC_RefGene_Name), ";", fixed = TRUE)
is_ox <- vapply(split_genes, function(z) any(z %in% ox_genes), logical(1))
is_promoter <- grepl("TSS1500|TSS200|1stExon",
                     as.character(ann$UCSC_RefGene_Group))
ox_promoter_probes <- ann$Probe[is_ox & is_promoter]
ox_promoter_probes <- intersect(ox_promoter_probes, rownames(M_cc))
writeLines(ox_promoter_probes,
           file.path(out_dir, "GSE298549_OXPHOS_promoter_probes.txt"))

# Competitive methylation gene-set testing with probe-number bias correction.
# Multiple cutoffs are sensitivity analyses; no cutoff is selected by its result.
collection <- list(GO_0006119_OXPHOS = ox_genes)
# gsameth expects Entrez IDs, rather than gene symbols.
ox_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = ox_genes, keytype = "SYMBOL", columns = "ENTREZID"))
ox_entrez <- unique(na.omit(as.character(ox_map$ENTREZID)))
collection <- list(GO_0006119_OXPHOS = ox_entrez)
gs_res <- rbindlist(lapply(c(0.001, 0.005, 0.01), function(cutoff) {
  sig <- ewas$Probe[ewas$P.Value < cutoff]
  if (length(sig) < 2) {
    return(data.table(CpG_P_cutoff = cutoff, N_selected_CpGs = length(sig),
                      Status = "Too few CpGs for gsameth"))
  }
  z <- tryCatch(
    missMethyl::gsameth(sig.cpg = sig, all.cpg = rownames(M_cc),
                        collection = collection, array.type = "450K",
                        plot.bias = FALSE),
    error = function(e) e)
  if (inherits(z, "error")) {
    return(data.table(CpG_P_cutoff = cutoff, N_selected_CpGs = length(sig),
                      Status = conditionMessage(z)))
  }
  ans <- as.data.table(z, keep.rownames = "GeneSet")
  ans[, `:=`(CpG_P_cutoff = cutoff, N_selected_CpGs = length(sig), Status = "OK")]
  ans
}), fill = TRUE)
fwrite(gs_res, file.path(out_dir, "GSE298549_OXPHOS_gsameth.tsv"), sep = "\t")

# A transparent sample-level targeted sensitivity score: standardize each
# OXPHOS promoter CpG across people, then average.  Direction is global promoter
# methylation only and must not be interpreted as gene-expression direction.
if (length(ox_promoter_probes) >= 10) {
  ox_m <- M_cc[ox_promoter_probes, , drop = FALSE]
  ox_z <- t(scale(t(ox_m)))
  score <- colMeans(ox_z, na.rm = TRUE)
  score <- as.numeric(scale(score))
  score_dat <- cbind(meta_cc, OXPHOS_promoter_score_z = score)
  score_fit <- lm(OXPHOS_promoter_score_z ~ WLE_z + Age_z + Sex + Position,
                  data = score_dat)
  score_ct <- lmtest::coeftest(score_fit,
                              vcov. = sandwich::vcovHC(score_fit, type = "HC3"))
  b <- unname(score_ct["WLE_z", 1]); se <- unname(score_ct["WLE_z", 2])
  score_result <- data.table(
    Analysis = "OXPHOS promoter mean-z methylation score",
    N = nobs(score_fit), N_CpGs = length(ox_promoter_probes),
    Beta_per_WLE_SD = b, SE_HC3 = se,
    CI_low = b - 1.96 * se, CI_high = b + 1.96 * se,
    P = unname(score_ct["WLE_z", 4]),
    Caveat = "Processed-data sensitivity; adjusted for age, sex, and array position; blood cell proportions unavailable")
  fwrite(score_result, file.path(out_dir, "GSE298549_OXPHOS_score_model.tsv"), sep = "\t")
  fwrite(as.data.table(score_dat)[, .(GSM, ArrayID, Age, Sex, WLE,
                                      OXPHOS_promoter_score_z)],
         file.path(out_dir, "GSE298549_OXPHOS_sample_scores.tsv"), sep = "\t")

  p_score <- ggplot(score_dat, aes(WLE, OXPHOS_promoter_score_z, color = Sex)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black") +
    theme_bw(base_size = 12) +
    labs(title = "GSE298549 MDD: life events and OXPHOS promoter methylation",
         subtitle = "Unadjusted visualization; inference uses the covariate-adjusted HC3 model",
         x = "Weighted life-events score", y = "Standardized OXPHOS promoter methylation score")
  ggsave(file.path(fig_dir, "GSE298549_WLE_OXPHOS_score.png"), p_score,
         width = 8.5, height = 6.2, dpi = 300)
  ggsave(file.path(fig_dir, "GSE298549_WLE_OXPHOS_score.pdf"), p_score,
         width = 8.5, height = 6.2)
}

# Genome-wide diagnostic plots.
plot_tab <- as.data.table(ewas)
plot_tab[, minus_log10_p := -log10(pmax(P.Value, .Machine$double.xmin))]
plot_tab[, Label := ifelse(Probe %in% reported, Probe, NA_character_)]
p_volcano <- ggplot(plot_tab, aes(logFC, minus_log10_p)) +
  geom_point(color = "grey65", alpha = 0.35, size = 0.7) +
  geom_point(data = plot_tab[Probe %in% reported], color = "#C0392B", size = 2) +
  geom_text_repel(aes(label = Label), na.rm = TRUE, max.overlaps = Inf) +
  theme_bw(base_size = 12) +
  labs(title = "GSE298549 weighted life-events EWAS in MDD",
       subtitle = "Effect is M-value change per 1 SD higher weighted life-events score",
       x = "Adjusted M-value change", y = expression(-log[10](P)))
ggsave(file.path(fig_dir, "GSE298549_wLE_volcano.png"), p_volcano,
       width = 8.5, height = 6.4, dpi = 300)
ggsave(file.path(fig_dir, "GSE298549_wLE_volcano.pdf"), p_volcano,
       width = 8.5, height = 6.4)

obs <- sort(ewas$P.Value)
exp_p <- ppoints(length(obs))
lambda_gc <- median(qchisq(1 - obs, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
qq <- data.table(Expected = -log10(exp_p), Observed = -log10(obs))
p_qq <- ggplot(qq, aes(Expected, Observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
  geom_point(alpha = 0.35, size = 0.8) +
  theme_bw(base_size = 12) +
  labs(title = "GSE298549 weighted life-events EWAS Q-Q plot",
       subtitle = sprintf("Genomic inflation factor lambda = %.3f", lambda_gc),
       x = expression(Expected~~-log[10](P)), y = expression(Observed~~-log[10](P)))
ggsave(file.path(fig_dir, "GSE298549_wLE_QQ.png"), p_qq,
       width = 6.6, height = 6.2, dpi = 300)
ggsave(file.path(fig_dir, "GSE298549_wLE_QQ.pdf"), p_qq,
       width = 6.6, height = 6.2)

# Optional DMR step for machines with adequate RAM.  Raw-IDAT analysis with cell
# adjustment remains preferable for final inference.
if (tolower(Sys.getenv("RUN_DMR", "false")) %in% c("true", "1", "yes")) {
  if (!requireNamespace("DMRcate", quietly = TRUE)) stop("DMRcate is not installed")
  dmr_ann <- DMRcate::cpg.annotate(
    datatype = "array", object = M_cc, what = "M",
    analysis.type = "differential", design = design, coef = "WLE_z",
    arraytype = "450K")
  dmr <- DMRcate::dmrcate(dmr_ann, lambda = 1000, C = 2)
  dmr_ranges <- DMRcate::extractRanges(dmr, genome = "hg19")
  fwrite(as.data.table(as.data.frame(dmr_ranges)),
         file.path(out_dir, "GSE298549_wLE_DMRcate.tsv"), sep = "\t")
}

summary_lines <- c(
  "GSE298549 processed-data sensitivity analysis",
  "=============================================",
  sprintf("MDD samples in processed matrix: %d", ncol(M)),
  sprintf("Complete weighted-life-events model: %d", ncol(M_cc)),
  sprintf("Autosomal probes tested: %d", nrow(M_cc)),
  sprintf("FDR-significant CpGs: %d", sum(ewas$adj.P.Val < 0.05, na.rm = TRUE)),
  sprintf("Lambda GC: %.4f", lambda_gc),
  "Primary coefficient: M-value change per 1 SD higher weighted life-events score",
  "Primary covariates: age, sex, array position",
  sprintf("Sparse-slide sensitivity design coefficients: %d", ncol(design_slide)),
  "Limitation: public processed-data analysis lacks estimated blood cell proportions.")
writeLines(summary_lines, file.path(out_dir, "GSE298549_analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "GSE298549_sessionInfo.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

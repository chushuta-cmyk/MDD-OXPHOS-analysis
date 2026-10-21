#!/usr/bin/env Rscript

# Independent anxiety-dimension validation of the pre-defined OXPHOS signal.
# Run from the MDD-MultiOmics repository root:
#   Rscript analysis_framework/12_GSE61672_anxiety_OXPHOS.R

repo_root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo_root, "analysis_framework", "_libpaths.R"))

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(AnnotationDbi)
  library(illuminaHumanv4.db)
  library(org.Hs.eg.db)
  library(limma)
  library(fgsea)
  library(singscore)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(sandwich)
  library(lmtest)
})

set.seed(61672)

base_dir <- file.path(repo_root, "analysis_framework", "GSE61672")
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "figures")
log_dir <- file.path(base_dir, "logs")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(data_dir, "GSE61672_series_matrix.txt.gz")
if (!file.exists(matrix_file)) {
  stop("Missing series matrix: ", matrix_file,
       "\nDownload it with analysis_framework/download_next_datasets.sh")
}

cat("===== GSE61672 anxiety/OXPHOS validation =====\n")
cat("Loading cached GEO series matrix\n")
gse <- getGEO(filename = matrix_file, getGPL = FALSE)
expr_probe <- exprs(gse)
pd <- pData(gse)
stopifnot(identical(colnames(expr_probe), rownames(pd)))

# GEO characteristics are position-shifted because not every sample has an
# anxiety case/control field. Parse key-value pairs across all characteristics.
char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
char_mat <- as.matrix(pd[, char_cols, drop = FALSE])
extract_key <- function(key) {
  prefix <- paste0(tolower(key), ":")
  apply(char_mat, 1, function(z) {
    z <- as.character(z)
    hit <- which(startsWith(tolower(trimws(z)), prefix))
    if (!length(hit)) return(NA_character_)
    trimws(sub("^[^:]+:", "", z[hit[1]]))
  })
}

meta <- data.frame(
  Sample = rownames(pd),
  Title = pd$title,
  Age = suppressWarnings(as.numeric(extract_key("age"))),
  Sex = extract_key("sex"),
  BMI = suppressWarnings(as.numeric(extract_key("body mass index"))),
  Ethnicity = extract_key("ethnicity"),
  GAD7 = suppressWarnings(as.numeric(extract_key("gad7 score"))),
  AnxietyGroupRaw = tolower(extract_key("anxiety case/control")),
  Batch = extract_key("hybridization batch"),
  RIN = suppressWarnings(as.numeric(extract_key("rin"))),
  stringsAsFactors = FALSE
)
meta$Sex <- factor(toupper(meta$Sex))
meta$Ethnicity3 <- factor(ifelse(meta$Ethnicity == "CAU", "CAU",
                                ifelse(meta$Ethnicity == "AFR", "AFR", "Other")))
meta$Batch <- factor(meta$Batch)
meta$AnxietyGroup <- factor(ifelse(meta$AnxietyGroupRaw == "case", "Case",
                                  ifelse(meta$AnxietyGroupRaw == "control", "Control", NA)),
                            levels = c("Control", "Case"))
rownames(meta) <- meta$Sample
write.csv(meta, file.path(data_dir, "GSE61672_metadata_parsed.csv"), row.names = FALSE)

cat("Expression:", nrow(expr_probe), "probes x", ncol(expr_probe), "samples\n")
cat("Valid numeric GAD7:", sum(!is.na(meta$GAD7)), "\n")
cat("Anxiety case/control subset:\n")
print(table(meta$AnxietyGroup, useNA = "ifany"))
cat("Parsed sex:\n")
print(table(meta$Sex, useNA = "ifany"))

# Probe-to-gene annotation and deterministic collapse: retain the most highly
# expressed probe for each symbol, matching the strategy used in prior cohorts.
symbol <- mapIds(illuminaHumanv4.db, keys = rownames(expr_probe),
                 keytype = "PROBEID", column = "SYMBOL", multiVals = "first")
ann <- data.table(Probe = rownames(expr_probe), Symbol = unname(symbol),
                  MeanExpr = rowMeans(expr_probe, na.rm = TRUE))
ann <- ann[!is.na(Symbol) & Symbol != ""]
setorder(ann, Symbol, -MeanExpr, Probe)
ann <- ann[!duplicated(Symbol)]
expr_gene <- expr_probe[ann$Probe, , drop = FALSE]
rownames(expr_gene) <- ann$Symbol
stopifnot(!anyDuplicated(rownames(expr_gene)))
saveRDS(expr_gene, file.path(data_dir, "GSE61672_gene_expression.rds"))
saveRDS(meta, file.path(data_dir, "GSE61672_metadata.rds"))

# Pre-defined GO:0006119 gene set; no outcome-driven selection in this cohort.
go_tab <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = "GO:0006119", keytype = "GOALL",
  columns = c("SYMBOL", "ONTOLOGYALL")))
ox_genes <- sort(unique(go_tab$SYMBOL[go_tab$ONTOLOGYALL == "BP" &
                                       !is.na(go_tab$SYMBOL)]))
ox_present <- intersect(ox_genes, rownames(expr_gene))
ox_nuclear <- ox_present[!grepl("^MT-", ox_present)]
core13 <- intersect(c("CHCHD2", "COX4I1", "COX5A", "COX5B", "COX6C",
                      "NDUFA13", "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFC2",
                      "NDUFS3", "SIRT3", "UQCR11"), rownames(expr_gene))
writeLines(ox_present, file.path(out_dir, "GSE61672_GO0006119_genes_present.txt"))

ranked <- rankGenes(expr_gene)
meta$Score_OXPHOS <- simpleScore(ranked, upSet = ox_present,
                                 centerScore = TRUE)$TotalScore
meta$Score_nuclear <- simpleScore(ranked, upSet = ox_nuclear,
                                  centerScore = TRUE)$TotalScore
meta$Score_core13 <- simpleScore(ranked, upSet = core13,
                                 centerScore = TRUE)$TotalScore
meta$Score_OXPHOS_z <- as.numeric(scale(meta$Score_OXPHOS))
meta$Score_nuclear_z <- as.numeric(scale(meta$Score_nuclear))
meta$Score_core13_z <- as.numeric(scale(meta$Score_core13))

znum <- function(x) as.numeric(scale(x))
meta$GAD7_z <- znum(meta$GAD7)
meta$Age_z <- znum(meta$Age)
meta$BMI_z <- znum(meta$BMI)
meta$RIN_z <- znum(meta$RIN)

required_primary <- c("GAD7_z", "Age_z", "BMI_z", "Ethnicity3", "Batch", "RIN_z")
required_sex <- c(required_primary, "Sex")
valid_cont <- complete.cases(meta[, required_primary])
valid_cont_sex <- complete.cases(meta[, required_sex])
valid_cc <- valid_cont & !is.na(meta$AnxietyGroup)
valid_cc_sex <- valid_cont_sex & !is.na(meta$AnxietyGroup)
cat("Complete continuous primary model N:", sum(valid_cont), "\n")
cat("Complete continuous sex-adjusted sensitivity N:", sum(valid_cont_sex), "\n")
cat("Complete case/control primary model N:", sum(valid_cc), "\n")
cat("Complete case/control sex-adjusted sensitivity N:", sum(valid_cc_sex), "\n")

hc3_term <- function(formula, data, term, analysis, outcome) {
  data <- droplevels(data)
  fit <- lm(formula, data = data)
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  if (!term %in% rownames(ct)) stop("Term not found: ", term)
  est <- unname(ct[term, 1]); se <- unname(ct[term, 2]); p <- unname(ct[term, 4])
  data.table(Analysis = analysis, Outcome = outcome, Term = term,
             N = nobs(fit), Beta = est, SE_HC3 = se,
             CI_low = est - 1.96 * se, CI_high = est + 1.96 * se,
             P = p, Model = paste(deparse(formula), collapse = ""))
}

cont_dat <- meta[valid_cont, , drop = FALSE]
cont_sex_dat <- meta[valid_cont_sex, , drop = FALSE]
cc_dat <- meta[valid_cc, , drop = FALSE]
cc_sex_dat <- meta[valid_cc_sex, , drop = FALSE]
score_results <- rbindlist(list(
  hc3_term(Score_OXPHOS_z ~ GAD7_z + Age_z + BMI_z + Ethnicity3 + Batch + RIN_z,
           cont_dat, "GAD7_z", "Continuous GAD-7 primary", "GO:0006119"),
  hc3_term(Score_nuclear_z ~ GAD7_z + Age_z + BMI_z + Ethnicity3 + Batch + RIN_z,
           cont_dat, "GAD7_z", "Continuous GAD-7 nuclear-only", "Nuclear OXPHOS"),
  hc3_term(Score_core13_z ~ GAD7_z + Age_z + BMI_z + Ethnicity3 + Batch + RIN_z,
           cont_dat, "GAD7_z", "Continuous GAD-7 core13", "Post-hoc core13"),
  hc3_term(Score_OXPHOS_z ~ AnxietyGroup + Age_z + BMI_z + Ethnicity3 + Batch + RIN_z,
           cc_dat, "AnxietyGroupCase", "Anxiety case-control primary", "GO:0006119"),
  hc3_term(Score_nuclear_z ~ AnxietyGroup + Age_z + BMI_z + Ethnicity3 + Batch + RIN_z,
           cc_dat, "AnxietyGroupCase", "Anxiety case-control nuclear-only", "Nuclear OXPHOS"),
  hc3_term(Score_core13_z ~ AnxietyGroup + Age_z + BMI_z + Ethnicity3 + Batch + RIN_z,
           cc_dat, "AnxietyGroupCase", "Anxiety case-control core13", "Post-hoc core13"),
  hc3_term(Score_OXPHOS_z ~ GAD7_z + Age_z + Sex + BMI_z + Ethnicity3 + Batch + RIN_z,
           cont_sex_dat, "GAD7_z", "Continuous GAD-7 sex-adjusted sensitivity", "GO:0006119"),
  hc3_term(Score_OXPHOS_z ~ AnxietyGroup + Age_z + Sex + BMI_z + Ethnicity3 + Batch + RIN_z,
           cc_sex_dat, "AnxietyGroupCase", "Anxiety case-control sex-adjusted sensitivity", "GO:0006119")
))

# Sex-interaction sensitivity, explicitly tested rather than comparing stratum P values.
score_results <- rbind(score_results,
  hc3_term(Score_OXPHOS_z ~ GAD7_z * Sex + Age_z + BMI_z + Ethnicity3 + Batch + RIN_z,
           cont_sex_dat, "GAD7_z:SexM", "Continuous GAD-7 sex interaction", "GO:0006119"),
  fill = TRUE)
fwrite(score_results, file.path(out_dir, "GSE61672_OXPHOS_score_models.tsv"), sep = "\t")

run_limma <- function(valid, design_formula, coef_name, label) {
  d <- droplevels(meta[valid, , drop = FALSE])
  design <- model.matrix(design_formula, data = d)
  stopifnot(coef_name %in% colnames(design))
  if (qr(design)$rank < ncol(design)) {
    nonest <- limma::nonEstimable(design)
    if (coef_name %in% nonest) stop("Primary coefficient is non-estimable: ", coef_name)
    design <- design[, !colnames(design) %in% nonest, drop = FALSE]
  }
  y <- expr_gene[, rownames(d), drop = FALSE]
  fit <- eBayes(lmFit(y, design), robust = TRUE, trend = TRUE)
  tab <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  tab$GeneSymbol <- rownames(tab)
  tab <- tab[, c("GeneSymbol", setdiff(colnames(tab), "GeneSymbol"))]
  fwrite(as.data.table(tab), file.path(out_dir, paste0("GSE61672_DE_", label, ".tsv")), sep = "\t")

  idx <- match(ox_present, rownames(y)); idx <- idx[!is.na(idx)]
  cam <- camera(y, index = list(GO_0006119_OXPHOS = idx), design = design,
                contrast = which(colnames(design) == coef_name))
  cam$GeneSet <- rownames(cam)
  cam$Analysis <- label

  stats <- fit$t[, coef_name]
  names(stats) <- rownames(fit$t)
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  fg <- fgseaMultilevel(pathways = list(GO_0006119_OXPHOS = ox_present),
                       stats = stats, minSize = 10, maxSize = 500,
                       eps = 0)
  fg$Analysis <- label
  list(table = tab, camera = as.data.table(cam), fgsea = as.data.table(fg),
       fit = fit, design = design)
}

lim_cont <- run_limma(valid_cont,
  ~ GAD7_z + Age_z + BMI_z + Ethnicity3 + Batch + RIN_z,
  "GAD7_z", "continuous_GAD7")
lim_cc <- run_limma(valid_cc,
  ~ AnxietyGroup + Age_z + BMI_z + Ethnicity3 + Batch + RIN_z,
  "AnxietyGroupCase", "case_vs_control")
fwrite(rbindlist(list(lim_cont$camera, lim_cc$camera), fill = TRUE),
       file.path(out_dir, "GSE61672_OXPHOS_CAMERA.tsv"), sep = "\t")
fwrite(rbindlist(list(lim_cont$fgsea, lim_cc$fgsea), fill = TRUE),
       file.path(out_dir, "GSE61672_OXPHOS_fgsea.tsv"), sep = "\t")

# Sample-level table for audit and later meta-analysis.
score_cols <- c("Sample", "Age", "Sex", "BMI", "Ethnicity", "GAD7",
                "AnxietyGroup", "Batch", "RIN", "Score_OXPHOS",
                "Score_OXPHOS_z", "Score_nuclear_z", "Score_core13_z")
fwrite(as.data.table(meta[, score_cols]),
       file.path(out_dir, "GSE61672_sample_scores.tsv"), sep = "\t")

theme_set(theme_bw(base_size = 11))

p_scatter <- ggplot(cont_dat, aes(GAD7, Score_OXPHOS_z, colour = Sex)) +
  geom_point(alpha = 0.55, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black",
              linewidth = 0.8) +
  scale_colour_manual(values = c(F = "#CC6677", M = "#4477AA")) +
  labs(title = "GSE61672: anxiety severity and OXPHOS score",
       subtitle = "Pre-defined GO:0006119 singscore; line is unadjusted visualization",
       x = "GAD-7 score", y = "Standardized OXPHOS score", colour = "Sex")
ggsave(file.path(fig_dir, "GSE61672_GAD7_OXPHOS_scatter.png"), p_scatter,
       width = 7.2, height = 5.2, dpi = 400)
ggsave(file.path(fig_dir, "GSE61672_GAD7_OXPHOS_scatter.pdf"), p_scatter,
       width = 7.2, height = 5.2)

p_box <- ggplot(cc_dat, aes(AnxietyGroup, Score_OXPHOS_z, fill = AnxietyGroup)) +
  geom_violin(trim = FALSE, alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.28, size = 1) +
  scale_fill_manual(values = c(Control = "#66C2A5", Case = "#FC8D62")) +
  guides(fill = "none") +
  labs(title = "GSE61672: anxiety case-control OXPHOS score",
       subtitle = "Plot is unadjusted; inference uses the covariate-adjusted HC3 model",
       x = NULL, y = "Standardized OXPHOS score")
ggsave(file.path(fig_dir, "GSE61672_anxiety_case_control_OXPHOS.png"), p_box,
       width = 6.5, height = 5.2, dpi = 400)
ggsave(file.path(fig_dir, "GSE61672_anxiety_case_control_OXPHOS.pdf"), p_box,
       width = 6.5, height = 5.2)

forest <- score_results[grepl("primary|nuclear", Analysis, ignore.case = TRUE)]
forest[, Analysis := factor(Analysis, levels = rev(Analysis))]
p_forest <- ggplot(forest, aes(Beta, Analysis)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.18,
                 colour = "#1B4F72") +
  geom_point(size = 2.6, colour = "#C0392B") +
  labs(title = "GSE61672 OXPHOS sensitivity models",
       subtitle = "Adjusted standardized effect with HC3 95% confidence interval",
       x = "Adjusted standardized effect", y = NULL)
ggsave(file.path(fig_dir, "GSE61672_OXPHOS_sensitivity_forest.png"), p_forest,
       width = 8.2, height = 4.8, dpi = 400)
ggsave(file.path(fig_dir, "GSE61672_OXPHOS_sensitivity_forest.pdf"), p_forest,
       width = 8.2, height = 4.8)

vol <- as.data.table(lim_cont$table)
vol[, neglog10P := -log10(pmax(P.Value, .Machine$double.xmin))]
vol[, IsOXPHOS := GeneSymbol %in% ox_present]
lab <- vol[IsOXPHOS == TRUE][order(P.Value)][1:min(12, .N)]
p_vol <- ggplot(vol, aes(logFC, neglog10P)) +
  geom_point(aes(colour = IsOXPHOS), alpha = 0.55, size = 1) +
  scale_colour_manual(values = c(`FALSE` = "grey72", `TRUE` = "#D55E00")) +
  geom_text_repel(data = lab, aes(label = GeneSymbol), size = 3,
                  max.overlaps = Inf, box.padding = 0.3) +
  guides(colour = guide_legend(title = "GO:0006119")) +
  labs(title = "GSE61672 continuous GAD-7 transcriptome",
       subtitle = "Covariate-adjusted limma; labelled genes are pre-defined OXPHOS members",
       x = "Expression change per 1 SD higher GAD-7", y = "-log10(P)")
ggsave(file.path(fig_dir, "GSE61672_continuous_GAD7_volcano.png"), p_vol,
       width = 7.5, height = 5.7, dpi = 400)
ggsave(file.path(fig_dir, "GSE61672_continuous_GAD7_volcano.pdf"), p_vol,
       width = 7.5, height = 5.7)

# PCA audit using the 1,000 most variable genes.
gene_var <- apply(expr_gene, 1, var, na.rm = TRUE)
top_var <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(1000, length(gene_var)))]
pc <- prcomp(t(expr_gene[top_var, , drop = FALSE]), center = TRUE, scale. = FALSE)
pcd <- data.frame(Sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                  GAD7 = meta[rownames(pc$x), "GAD7"],
                  Batch = meta[rownames(pc$x), "Batch"])
pcd$GAD7_plot <- ifelse(is.na(pcd$GAD7), 0, pcd$GAD7)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
p_pca <- ggplot(pcd, aes(PC1, PC2, colour = Batch, size = GAD7_plot)) +
  geom_point(alpha = 0.65) +
  scale_size_continuous(range = c(1, 4), name = "GAD-7") +
  labs(title = "GSE61672 expression PCA",
       x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]))
ggsave(file.path(fig_dir, "GSE61672_expression_PCA.png"), p_pca,
       width = 7.2, height = 5.5, dpi = 400)

summary_lines <- c(
  "GSE61672 independent anxiety/OXPHOS validation",
  sprintf("Samples in GEO matrix: %d", ncol(expr_probe)),
  sprintf("Genes after probe collapse: %d", nrow(expr_gene)),
  sprintf("GO:0006119 genes present: %d", length(ox_present)),
  sprintf("Continuous GAD-7 primary complete-case N: %d", sum(valid_cont)),
  sprintf("Continuous GAD-7 sex-adjusted sensitivity N: %d", sum(valid_cont_sex)),
  sprintf("Case/control primary complete-case N: %d", sum(valid_cc)),
  sprintf("Case/control sex-adjusted sensitivity N: %d", sum(valid_cc_sex)),
  "",
  capture.output(print(score_results)),
  "",
  "CAMERA:", capture.output(print(rbindlist(list(lim_cont$camera, lim_cc$camera), fill = TRUE))),
  "",
  "fgsea:", capture.output(print(rbindlist(list(lim_cont$fgsea, lim_cc$fgsea), fill = TRUE))),
  "",
  "Interpretation guardrail: GSE61672 tests an anxiety dimension in a community blood cohort;",
  "it is not a direct independent replication of the MDD+GAD versus MDD-only contrast."
)
writeLines(summary_lines, file.path(out_dir, "GSE61672_analysis_summary.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("===== GSE61672 pipeline complete =====\n")
cat("Results:", out_dir, "\nFigures:", fig_dir, "\n")

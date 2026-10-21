# ==============================================================================
# 11_blood_OXPHOS_strict_sensitivity_meta.R
#
# Public-data sensitivity analysis and cross-cohort integration of a fixed
# oxidative-phosphorylation gene set (GO:0006119) in four MDD blood datasets.
#
# Primary evidence set: GSE98793, GSE39653, GSE19738.
# Sensitivity-only cohort: GSE251778, because public GEO metadata omit the age,
# BMI, RIN and blood-cell covariates used by the source publication and sex is
# perfectly confounded with sequencing platform/library protocol.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(singscore)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_arg)) dirname(normalizePath(script_arg)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."))
out_dir <- file.path(script_dir, "results_strict")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260920)

go_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = "GO:0006119", keytype = "GOALL", columns = "SYMBOL"))
oxphos_all <- sort(unique(na.omit(toupper(go_map$SYMBOL))))
oxphos_nuclear <- oxphos_all[!grepl("^MT-", oxphos_all)]
core13 <- c("CHCHD2", "COX4I1", "COX5A", "COX5B", "COX6C", "NDUFA13",
            "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFC2", "NDUFS3", "SIRT3",
            "UQCR11")

collapse_expression <- function(x) {
  stopifnot(is.matrix(x), !is.null(rownames(x)))
  genes <- toupper(trimws(rownames(x)))
  keep <- !is.na(genes) & nzchar(genes) & rowSums(is.finite(x)) == ncol(x)
  x <- x[keep, , drop = FALSE]
  genes <- genes[keep]
  if (anyDuplicated(genes)) {
    rv <- matrixStats::rowVars(x)
    ord <- order(genes, -rv)
    x <- x[ord, , drop = FALSE]
    genes <- genes[ord]
    first <- !duplicated(genes)
    x <- x[first, , drop = FALSE]
    genes <- genes[first]
  }
  rownames(x) <- genes
  x
}

score_expression <- function(expr, genes) {
  expr <- collapse_expression(expr)
  use <- intersect(genes, rownames(expr))
  if (length(use) < 10) stop("Fewer than 10 OXPHOS genes were measured.")
  ranked <- singscore::rankGenes(expr)
  sc <- singscore::simpleScore(ranked, upSet = use, centerScore = TRUE,
                               knownDirection = TRUE)$TotalScore
  sc <- as.numeric(sc)
  names(sc) <- colnames(expr)
  list(raw = sc, z = as.numeric(scale(sc)), measured = length(use))
}

hc3_term <- function(formula, data, term, cohort, analysis, gene_set,
                     caveat = "") {
  fit <- lm(formula, data = data)
  X <- model.matrix(fit)
  if (!term %in% colnames(X)) stop("Term not found: ", term)
  e <- residuals(fit)
  h <- hatvalues(fit)
  omega <- (e / pmax(1 - h, 1e-8))^2
  bread <- solve(crossprod(X))
  meat <- crossprod(X, X * omega)
  vc <- bread %*% meat %*% bread
  j <- match(term, colnames(X))
  beta <- coef(fit)[j]
  se <- sqrt(vc[j, j])
  z <- beta / se
  data.table(Cohort = cohort, Analysis = analysis, GeneSet = gene_set,
             N = nobs(fit), Beta = beta, SE = se,
             CI_low = beta - 1.96 * se, CI_high = beta + 1.96 * se,
             Z = z, P = 2 * pnorm(abs(z), lower.tail = FALSE),
             Term = term, Caveat = caveat)
}

prepare_score_data <- function(expr, meta, group, sex = NULL, age = NULL,
                               batch = NULL) {
  stopifnot(all(colnames(expr) %in% rownames(meta)))
  meta <- meta[colnames(expr), , drop = FALSE]
  s_all <- score_expression(expr, oxphos_all)
  s_nuc <- score_expression(expr, oxphos_nuclear)
  s_core <- score_expression(expr, core13)
  d <- data.frame(Sample = colnames(expr), Group = group(meta),
                  Score_all = s_all$z, Score_nuclear = s_nuc$z,
                  Score_core13 = s_core$z, stringsAsFactors = FALSE)
  d$Group <- factor(d$Group, levels = c("Control", "MDD"))
  if (!is.null(sex)) d$Sex <- factor(sex(meta), levels = c("Female", "Male"))
  if (!is.null(age)) d$Age <- as.numeric(age(meta))
  if (!is.null(batch)) d$Batch <- factor(batch(meta))
  attr(d, "measured") <- c(all = s_all$measured, nuclear = s_nuc$measured,
                            core13 = s_core$measured)
  d
}

score_results <- list()
sample_scores <- list()
camera_results <- list()

# ------------------------------------------------------------------------------
# GSE98793
# ------------------------------------------------------------------------------
expr98793 <- readRDS(file.path(repo_root,
  "transcriptomics/blood/GSE98793/results_fixed/GSE98793_expression_cleaned.rds"))
meta98793 <- readRDS(file.path(repo_root,
  "transcriptomics/blood/GSE98793/results_fixed/GSE98793_metadata.rds"))
meta98793 <- meta98793[match(colnames(expr98793), meta98793$Sample), , drop = FALSE]
rownames(meta98793) <- meta98793$Sample
d98793 <- prepare_score_data(
  expr98793, meta98793,
  group = function(m) m$disease_group,
  sex = function(m) ifelse(toupper(m$gender) == "F", "Female", "Male"),
  age = function(m) m$age,
  batch = function(m) m$batch)
d98793$GAD <- meta98793[d98793$Sample, "GAD"]
sample_scores$GSE98793 <- data.table(Cohort = "GSE98793", d98793)

score_results$g98793_main <- hc3_term(
  Score_all ~ Group + Age + Sex + Batch, d98793, "GroupMDD", "GSE98793",
  "Primary covariate-adjusted", "GO:0006119")
score_results$g98793_nuc <- hc3_term(
  Score_nuclear ~ Group + Age + Sex + Batch, d98793, "GroupMDD", "GSE98793",
  "Nuclear-encoded sensitivity", "GO:0006119 nuclear only")
score_results$g98793_core <- hc3_term(
  Score_core13 ~ Group + Age + Sex + Batch, d98793, "GroupMDD", "GSE98793",
  "Post-hoc core-gene score", "13-gene core",
  "Discovery-derived score; descriptive only")
nogad <- d98793[d98793$Group == "Control" | d98793$GAD == 0, , drop = FALSE]
score_results$g98793_nogad <- hc3_term(
  Score_all ~ Group + Age + Sex + Batch, nogad, "GroupMDD", "GSE98793",
  "MDD without GAD vs control", "GO:0006119")
score_results$g98793_interaction <- hc3_term(
  Score_all ~ Group * Sex + Age + Batch, d98793, "GroupMDD:SexMale", "GSE98793",
  "Diagnosis by sex interaction", "GO:0006119",
  "Positive interaction means a less negative or more positive MDD effect in males")
for (sx in levels(d98793$Sex)) {
  z <- droplevels(d98793[d98793$Sex == sx, ])
  score_results[[paste0("g98793_", sx)]] <- hc3_term(
    Score_all ~ Group + Age + Batch, z, "GroupMDD", "GSE98793",
    paste0(sx, " stratum"), "GO:0006119")
}
design98793 <- model.matrix(~ Group + Age + Sex + Batch, data = d98793)
cam98793 <- camera(collapse_expression(expr98793),
                   index = list(OXPHOS = intersect(oxphos_all,
                     rownames(collapse_expression(expr98793)))),
                   design = design98793, contrast = "GroupMDD")
camera_results$GSE98793 <- data.table(Cohort = "GSE98793", cam98793,
                                      Model = "Group + age + sex + batch")

# ------------------------------------------------------------------------------
# GSE39653
# ------------------------------------------------------------------------------
expr39653 <- readRDS(file.path(repo_root,
  "transcriptomics/blood/GSE39653/results_fixed/GSE39653_gene_expression.rds"))
meta39653 <- readRDS(file.path(repo_root,
  "transcriptomics/blood/GSE39653/results_fixed/GSE39653_metadata.rds"))
meta39653 <- meta39653[match(colnames(expr39653), rownames(meta39653)), , drop = FALSE]
rownames(meta39653) <- colnames(expr39653)
d39653 <- prepare_score_data(expr39653, meta39653,
                             group = function(m) as.character(m$Diagnosis))
sample_scores$GSE39653 <- data.table(Cohort = "GSE39653", d39653)
score_results$g39653_main <- hc3_term(
  Score_all ~ Group, d39653, "GroupMDD", "GSE39653",
  "Primary unadjusted", "GO:0006119",
  "No additional individual-level covariates available in the local GEO metadata")
score_results$g39653_nuc <- hc3_term(
  Score_nuclear ~ Group, d39653, "GroupMDD", "GSE39653",
  "Nuclear-encoded sensitivity", "GO:0006119 nuclear only")
score_results$g39653_core <- hc3_term(
  Score_core13 ~ Group, d39653, "GroupMDD", "GSE39653",
  "Post-hoc core-gene score", "13-gene core",
  "Discovery-derived score; descriptive only")
design39653 <- model.matrix(~ Group, data = d39653)
expr39653c <- collapse_expression(expr39653)
cam39653 <- camera(expr39653c,
                   index = list(OXPHOS = intersect(oxphos_all, rownames(expr39653c))),
                   design = design39653, contrast = "GroupMDD")
camera_results$GSE39653 <- data.table(Cohort = "GSE39653", cam39653,
                                      Model = "Group")

# ------------------------------------------------------------------------------
# GSE19738
# ------------------------------------------------------------------------------
expr19738 <- readRDS(file.path(repo_root,
  "analysis_framework/GSE19738/results/GSE19738_gene_expression.rds"))
meta19738 <- readRDS(file.path(repo_root,
  "analysis_framework/GSE19738/results/GSE19738_metadata.rds"))
meta19738 <- meta19738[colnames(expr19738), , drop = FALSE]
d19738 <- prepare_score_data(
  expr19738, meta19738,
  group = function(m) ifelse(as.character(m$MDD) == "1", "MDD", "Control"),
  sex = function(m) ifelse(tolower(as.character(m$Sex)) == "female", "Female", "Male"),
  age = function(m) m$Age)
sample_scores$GSE19738 <- data.table(Cohort = "GSE19738", d19738)
score_results$g19738_main <- hc3_term(
  Score_all ~ Group + Age + Sex, d19738, "GroupMDD", "GSE19738",
  "Primary covariate-adjusted", "GO:0006119")
score_results$g19738_nuc <- hc3_term(
  Score_nuclear ~ Group + Age + Sex, d19738, "GroupMDD", "GSE19738",
  "Nuclear-encoded sensitivity", "GO:0006119 nuclear only")
score_results$g19738_core <- hc3_term(
  Score_core13 ~ Group + Age + Sex, d19738, "GroupMDD", "GSE19738",
  "Post-hoc core-gene score", "13-gene core",
  "Discovery-derived score; descriptive only")
score_results$g19738_interaction <- hc3_term(
  Score_all ~ Group * Sex + Age, d19738, "GroupMDD:SexMale", "GSE19738",
  "Diagnosis by sex interaction", "GO:0006119",
  "Positive interaction means a less negative or more positive MDD effect in males")
for (sx in levels(d19738$Sex)) {
  z <- droplevels(d19738[d19738$Sex == sx, ])
  score_results[[paste0("g19738_", sx)]] <- hc3_term(
    Score_all ~ Group + Age, z, "GroupMDD", "GSE19738",
    paste0(sx, " stratum"), "GO:0006119")
}
design19738 <- model.matrix(~ Group + Age + Sex, data = d19738)
expr19738c <- collapse_expression(expr19738)
cam19738 <- camera(expr19738c,
                   index = list(OXPHOS = intersect(oxphos_all, rownames(expr19738c))),
                   design = design19738, contrast = "GroupMDD")
camera_results$GSE19738 <- data.table(Cohort = "GSE19738", cam19738,
                                      Model = "Group + age + sex")

# ------------------------------------------------------------------------------
# GSE251778: public-data sensitivity analysis
# ------------------------------------------------------------------------------
count_file <- file.path(repo_root,
  "transcriptomics/blood/GSE251778/data/GSE251778_Mokhtari_rna_counts.csv.gz")
meta_file <- file.path(repo_root,
  "transcriptomics/blood/GSE251778/data/GSE251778_metadata.csv")
count_dt <- fread(count_file, check.names = FALSE)
gene_col <- if ("gene_ID" %in% names(count_dt)) "gene_ID" else names(count_dt)[2]
ens <- sub("\\..*$", "", as.character(count_dt[[gene_col]]))
sample_cols <- setdiff(names(count_dt), names(count_dt)[1:2])
counts <- as.matrix(count_dt[, ..sample_cols])
storage.mode(counts) <- "integer"
rownames(counts) <- ens
colnames(counts) <- sub("_rna$", "", sample_cols)
if (anyDuplicated(rownames(counts))) counts <- rowsum(counts, rownames(counts), reorder = FALSE)
sym <- unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(counts),
  keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
keep_sym <- !is.na(sym) & nzchar(sym)
counts_symbol <- rowsum(counts[keep_sym, , drop = FALSE],
                        toupper(sym[keep_sym]), reorder = FALSE)
meta251778 <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
meta251778 <- meta251778[match(colnames(counts_symbol), meta251778$SampleID), , drop = FALSE]
rownames(meta251778) <- meta251778$SampleID
meta251778$Group <- factor(meta251778$Group, levels = c("Control", "MDD"))
meta251778$Sex <- factor(meta251778$Sex, levels = c("Female", "Male"))

gsea_one <- function(stats) {
  stats <- sort(stats[is.finite(stats) & !duplicated(names(stats))], decreasing = TRUE)
  y <- fgseaMultilevel(list(OXPHOS = oxphos_all), stats = stats,
                       minSize = 10, maxSize = 500, eps = 0)
  data.table(NES = y$NES, P = y$pval, Size = y$size,
             LeadingEdge = paste(y$leadingEdge[[1]], collapse = ";"))
}

# GSE98793 comorbid-anxiety sensitivity using the already fitted contrasts.
# The overall MDD result combines equal-sized MDD-with-GAD and MDD-without-GAD
# groups; separating them tests whether the pathway signal is diagnosis-wide.
gad_files <- c(
  MDD_without_GAD_vs_control = file.path(repo_root,
    "transcriptomics/blood/GSE98793/results_fixed/GSE98793_DEG_MDDnoGAD_vs_Ctrl.csv"),
  MDD_with_GAD_vs_control = file.path(repo_root,
    "transcriptomics/blood/GSE98793/results_fixed/GSE98793_DEG_MDDGAD_vs_Ctrl.csv"),
  MDD_with_GAD_vs_without_GAD = file.path(repo_root,
    "transcriptomics/blood/GSE98793/results_fixed/GSE98793_DEG_MDDGAD_vs_MDDnoGAD.csv"))
gad_gsea <- rbindlist(lapply(names(gad_files), function(nm) {
  x <- fread(gad_files[[nm]])
  ranks <- setNames(x$t, x[[1]])
  y <- gsea_one(ranks)
  data.table(Contrast = nm, NES = y$NES, P = y$P, Size = y$Size)
}))
fwrite(gad_gsea, file.path(out_dir, "GSE98793_GAD_OXPHOS_sensitivity.tsv"), sep = "\t")

threshold_results <- list()
dds10 <- NULL
vst10 <- NULL
for (thr in c(5, 10, 20)) {
  keep <- rowMeans(counts_symbol) > thr
  dds <- DESeqDataSetFromMatrix(counts_symbol[keep, , drop = FALSE],
                                meta251778, design = ~ Sex + Group)
  dds <- DESeq(dds, quiet = TRUE)
  rr <- results(dds, contrast = c("Group", "MDD", "Control"))
  gs <- gsea_one(setNames(rr$stat, rownames(rr)))
  threshold_results[[as.character(thr)]] <- data.table(
    Cohort = "GSE251778", Method = "DESeq2", Filter = paste0("mean count >", thr),
    Genes = sum(keep), NES = gs$NES, P = gs$P)
  if (thr == 10) {
    dds10 <- dds
    vst10 <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  }
  rm(dds, rr); gc()
}

# Independent RNA-seq method sensitivity: TMM + voom + robust limma.
y251 <- DGEList(counts_symbol)
keep_voom <- filterByExpr(y251, design = model.matrix(~ Sex + Group, meta251778))
y251 <- calcNormFactors(y251[keep_voom, , keep.lib.sizes = FALSE])
design251 <- model.matrix(~ Sex + Group, meta251778)
v251 <- voom(y251, design251, plot = FALSE)
fit251 <- eBayes(lmFit(v251, design251), trend = TRUE, robust = TRUE)
t251 <- fit251$t[, "GroupMDD"]
gs251_voom <- gsea_one(setNames(t251, rownames(v251$E)))
threshold_results$voom <- data.table(
  Cohort = "GSE251778", Method = "edgeR limma-voom", Filter = "filterByExpr",
  Genes = nrow(v251$E), NES = gs251_voom$NES, P = gs251_voom$P)

d251778 <- prepare_score_data(
  vst10, meta251778,
  group = function(m) as.character(m$Group),
  sex = function(m) as.character(m$Sex))
sample_scores$GSE251778 <- data.table(Cohort = "GSE251778", d251778)
score_results$g251778_main <- hc3_term(
  Score_all ~ Group + Sex, d251778, "GroupMDD", "GSE251778",
  "Public-data primary", "GO:0006119",
  "Age, BMI, RIN and blood-cell covariates unavailable; sex and platform are confounded")
score_results$g251778_nuc <- hc3_term(
  Score_nuclear ~ Group + Sex, d251778, "GroupMDD", "GSE251778",
  "Nuclear-encoded sensitivity", "GO:0006119 nuclear only",
  "Age, BMI, RIN and blood-cell covariates unavailable")
score_results$g251778_core <- hc3_term(
  Score_core13 ~ Group + Sex, d251778, "GroupMDD", "GSE251778",
  "Post-hoc core-gene score", "13-gene core",
  "Discovery-derived and missing published covariates; descriptive only")
score_results$g251778_interaction <- hc3_term(
  Score_all ~ Group * Sex, d251778, "GroupMDD:SexMale", "GSE251778",
  "Diagnosis by sex interaction", "GO:0006119",
  "Sex is perfectly confounded with platform; interaction cannot be separated from platform-specific effects")
for (sx in levels(d251778$Sex)) {
  z <- droplevels(d251778[d251778$Sex == sx, ])
  score_results[[paste0("g251778_", sx)]] <- hc3_term(
    Score_all ~ Group, z, "GroupMDD", "GSE251778",
    paste0(sx, " stratum"), "GO:0006119",
    "Sex stratum is also a platform/library-protocol stratum")
}

# SVA sensitivity. OXPHOS genes are excluded when estimating latent variables
# to reduce direct removal of the outcome score.
if (requireNamespace("sva", quietly = TRUE)) {
  non_ox <- vst10[!rownames(vst10) %in% oxphos_all, , drop = FALSE]
  mod <- model.matrix(~ Sex + Group, meta251778)
  mod0 <- model.matrix(~ Sex, meta251778)
  nsv <- tryCatch(sva::num.sv(non_ox, mod, method = "be"), error = function(e) 0L)
  nsv <- min(as.integer(nsv), 10L)
  if (nsv > 0) {
    svobj <- sva::sva(non_ox, mod, mod0, n.sv = nsv)
    for (j in seq_len(nsv)) d251778[[paste0("SV", j)]] <- svobj$sv[, j]
    rhs <- paste(c("Group", "Sex", paste0("SV", seq_len(nsv))), collapse = " + ")
    score_results$g251778_sva <- hc3_term(
      as.formula(paste("Score_all ~", rhs)), d251778, "GroupMDD", "GSE251778",
      paste0("SVA-adjusted public-data sensitivity with ", nsv, " SVs"),
      "GO:0006119",
      "Surrogate variables estimated from non-OXPHOS genes; exact published covariates remain unavailable")
  }
}

cam251 <- camera(v251$E,
                 index = list(OXPHOS = intersect(oxphos_all, rownames(v251$E))),
                 design = design251, contrast = "GroupMDD")
camera_results$GSE251778 <- data.table(Cohort = "GSE251778", cam251,
                                       Model = "Sex/platform + group; voom")

# ------------------------------------------------------------------------------
# Save within-cohort sensitivity results
# ------------------------------------------------------------------------------
score_tab <- rbindlist(score_results, fill = TRUE)
setorder(score_tab, Cohort, Analysis)
fwrite(score_tab, file.path(out_dir, "blood_OXPHOS_singscore_sensitivity.tsv"), sep = "\t")
fwrite(rbindlist(sample_scores, fill = TRUE),
       file.path(out_dir, "blood_OXPHOS_sample_scores.tsv"), sep = "\t")
camera_tab <- rbindlist(camera_results, fill = TRUE)
fwrite(camera_tab, file.path(out_dir, "blood_OXPHOS_CAMERA_sensitivity.tsv"), sep = "\t")
threshold_tab <- rbindlist(threshold_results)
fwrite(threshold_tab, file.path(out_dir, "GSE251778_method_filter_sensitivity.tsv"), sep = "\t")

# ------------------------------------------------------------------------------
# Random-effects meta-analysis with REML and Hartung-Knapp confidence intervals
# ------------------------------------------------------------------------------
reml_meta <- function(dat, label) {
  dat <- dat[is.finite(Beta) & is.finite(SE) & SE > 0]
  yi <- dat$Beta; vi <- dat$SE^2; k <- length(yi)
  if (k < 2) stop("At least two cohorts are required.")
  nll <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu)^2))
  }
  upper <- max(1, stats::var(yi) * 100, max(vi) * 100)
  tau2 <- optimize(nll, c(0, upper))$minimum
  if (nll(0) <= nll(tau2) + 1e-8) tau2 <- 0
  w <- 1 / (vi + tau2)
  mu <- sum(w * yi) / sum(w)
  se_wald <- sqrt(1 / sum(w))
  qstar <- sum(w * (yi - mu)^2) / (k - 1)
  se_hk <- sqrt(max(1, qstar) / sum(w))
  tcrit <- qt(0.975, df = k - 1)
  p_hk <- 2 * pt(abs(mu / se_hk), df = k - 1, lower.tail = FALSE)
  wf <- 1 / vi
  muf <- sum(wf * yi) / sum(wf)
  Q <- sum(wf * (yi - muf)^2)
  I2 <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  pred_crit <- if (k > 2) qt(0.975, df = k - 2) else qnorm(0.975)
  data.table(Analysis = label, K = k, Pooled_beta = mu,
             SE_HK = se_hk, CI_low_HK = mu - tcrit * se_hk,
             CI_high_HK = mu + tcrit * se_hk, P_HK = p_hk,
             Tau2_REML = tau2, Q = Q, Q_df = k - 1,
             Q_P = pchisq(Q, df = k - 1, lower.tail = FALSE), I2 = I2,
             Prediction_low = mu - pred_crit * sqrt(tau2 + se_wald^2),
             Prediction_high = mu + pred_crit * sqrt(tau2 + se_wald^2))
}

primary_names <- c("GSE98793", "GSE39653", "GSE19738")
main_scores <- score_tab[Analysis %in% c("Primary covariate-adjusted",
                                         "Primary unadjusted",
                                         "Public-data primary") &
                           GeneSet == "GO:0006119"]
primary_scores <- main_scores[Cohort %in% primary_names]
all_scores <- main_scores[Cohort %in% c(primary_names, "GSE251778")]
meta_tab <- rbindlist(list(
  reml_meta(primary_scores, "Primary three-cohort model"),
  reml_meta(all_scores, "Four-cohort sensitivity including GSE251778")))

loo <- rbindlist(lapply(primary_names, function(drop) {
  reml_meta(primary_scores[Cohort != drop], paste0("Primary leave out ", drop))
}))
loo_all <- rbindlist(lapply(all_scores$Cohort, function(drop) {
  reml_meta(all_scores[Cohort != drop], paste0("Four-cohort leave out ", drop))
}))
meta_all <- rbindlist(list(meta_tab, loo, loo_all), fill = TRUE)
fwrite(meta_all, file.path(out_dir, "blood_OXPHOS_random_effects_meta.tsv"), sep = "\t")

# Formal interaction synthesis. GSE251778 remains sensitivity-only.
interaction_scores <- score_tab[Analysis == "Diagnosis by sex interaction"]
interaction_meta <- rbindlist(list(
  reml_meta(interaction_scores[Cohort %in% c("GSE98793", "GSE19738")],
            "Interaction primary cohorts with sex data"),
  reml_meta(interaction_scores,
            "Interaction sensitivity including sex-platform-confounded GSE251778")))
fwrite(interaction_meta, file.path(out_dir, "blood_OXPHOS_sex_interaction_meta.tsv"), sep = "\t")

# Weighted signed-Z integration of the independent ranked-list GSEA results.
rank_file <- file.path(script_dir, "results", "predefined_OXPHOS_across_transcriptomes.tsv")
rank_tab <- fread(rank_file)
rank_ids <- c(GSE98793 = "GSE98793_blood", GSE39653 = "GSE39653_PBMC",
              GSE19738 = "GSE19738_blood", GSE251778 = "GSE251778_overall")
n_total <- c(GSE98793 = 192, GSE39653 = 45, GSE19738 = 67, GSE251778 = 169)
signed_z_combine <- function(ids, label) {
  x <- rank_tab[match(unname(rank_ids[ids]), Dataset)]
  z <- sign(x$NES) * qnorm(pmax(x$P, 1e-300) / 2, lower.tail = FALSE)
  w <- sqrt(n_total[ids])
  zz <- sum(w * z) / sqrt(sum(w^2))
  data.table(Analysis = label, K = length(ids), Signed_Z = zz,
             P_two_sided = 2 * pnorm(abs(zz), lower.tail = FALSE),
             Cohorts = paste(ids, collapse = ";"))
}
stouffer_tab <- rbindlist(list(
  signed_z_combine(primary_names, "Primary three-cohort ranked-list integration"),
  signed_z_combine(c(primary_names, "GSE251778"),
                   "Four-cohort ranked-list sensitivity"),
  rbindlist(lapply(primary_names, function(drop)
    signed_z_combine(setdiff(primary_names, drop), paste0("Leave out ", drop))))))
fwrite(stouffer_tab, file.path(out_dir, "blood_OXPHOS_weighted_Stouffer.tsv"), sep = "\t")

# ------------------------------------------------------------------------------
# Core-gene random-effects meta-analysis
# ------------------------------------------------------------------------------
gene_files <- list(
  GSE98793 = file.path(repo_root, "transcriptomics/blood/GSE98793/results_fixed/GSE98793_all_DEG_results.csv"),
  GSE39653 = file.path(repo_root, "transcriptomics/blood/GSE39653/results_fixed/GSE39653_all_DEG_results.csv"),
  GSE19738 = file.path(repo_root, "analysis_framework/GSE19738/results/deg_GSE19738.csv"),
  GSE251778 = file.path(repo_root, "transcriptomics/blood/GSE251778/results_fixed/GSE251778_DEG_overall.csv"))

read_gene_effect <- function(f, cohort) {
  x <- fread(f)
  gene <- if ("GeneSymbol" %in% names(x)) x$GeneSymbol else x$gene
  eff <- if ("log2FoldChange" %in% names(x)) x$log2FoldChange else x$logFC
  if ("lfcSE" %in% names(x)) {
    se <- x$lfcSE
  } else {
    tt <- if ("t" %in% names(x)) x$t else x$stat
    se <- abs(eff / tt)
  }
  data.table(Cohort = cohort, Gene = toupper(gene), Beta = eff, SE = se)[
    Gene %in% core13 & is.finite(Beta) & is.finite(SE) & SE > 0]
}
gene_effects <- rbindlist(Map(read_gene_effect, gene_files, names(gene_files)))
fwrite(gene_effects, file.path(out_dir, "core13_gene_effects_for_meta.tsv"), sep = "\t")
gene_meta <- rbindlist(lapply(core13, function(g) {
  a <- reml_meta(gene_effects[Gene == g & Cohort %in% primary_names],
                 "Primary three cohorts")
  b <- reml_meta(gene_effects[Gene == g], "Including GSE251778 sensitivity")
  a[, Gene := g]; b[, Gene := g]
  rbind(a, b, fill = TRUE)
}))
gene_meta[, FDR_HK := p.adjust(P_HK, method = "BH"), by = Analysis]
setcolorder(gene_meta, c("Gene", setdiff(names(gene_meta), "Gene")))
fwrite(gene_meta, file.path(out_dir, "core13_random_effects_meta.tsv"), sep = "\t")

# ------------------------------------------------------------------------------
# Main figures
# ------------------------------------------------------------------------------
forest_data <- copy(all_scores)
forest_data[, Label := Cohort]
pooled_plot <- meta_tab[, .(Cohort = Analysis, Analysis = Analysis,
  GeneSet = "GO:0006119", N = NA_integer_, Beta = Pooled_beta, SE = SE_HK,
  CI_low = CI_low_HK, CI_high = CI_high_HK, Label = Analysis)]
forest_plot <- rbindlist(list(forest_data, pooled_plot), fill = TRUE)
forest_plot$Label <- factor(forest_plot$Label, levels = rev(c(
  "GSE98793", "GSE39653", "GSE19738", "GSE251778",
  "Primary three-cohort model",
  "Four-cohort sensitivity including GSE251778")))
p_forest <- ggplot(forest_plot, aes(Beta, Label)) +
  geom_vline(xintercept = 0, color = "grey45", linetype = 2) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.18,
                 color = "#243B53") +
  geom_point(aes(shape = grepl("cohort", as.character(Label), ignore.case = TRUE)),
             size = 3, color = "#005A8D") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  labs(x = "Adjusted difference in standardized single-sample OXPHOS score\nMDD minus control",
       y = NULL, title = "Peripheral blood OXPHOS score effects",
       subtitle = "Cohort estimates use HC3 standard errors\nPooled estimates use REML and Hartung Knapp intervals") +
  theme_minimal(base_size = 11) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "Figure_blood_OXPHOS_score_forest.png"), p_forest,
       width = 9.2, height = 5.6, dpi = 320)

rank_plot <- copy(rank_tab)
rank_plot[, Tissue := ifelse(grepl("blood|PBMC|GSE251778|GSE19738", Dataset),
                             "Peripheral blood", "Brain")]
rank_plot[, Significant := P < 0.05]
p_rank <- ggplot(rank_plot, aes(NES, reorder(Dataset, NES), color = Significant)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 2.7) + facet_grid(Tissue ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey55")) +
  labs(x = "Normalized enrichment score for GO:0006119", y = NULL,
       color = "Nominal P < 0.05", title = "OXPHOS enrichment is consistent in three blood cohorts but heterogeneous across tissues") +
  theme_minimal(base_size = 10) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "Figure_OXPHOS_NES_all_transcriptomes.png"), p_rank,
       width = 9.2, height = 7.8, dpi = 320)

core_detail <- fread(file.path(script_dir, "results", "blood_OXPHOS_core_gene_cross_dataset.tsv"))
heat_ids <- c("GSE98793_blood", "GSE39653_PBMC", "GSE19738_blood",
              "GSE19738_female", "GSE19738_male",
              "GSE251778_female", "GSE251778_male")
heat <- dcast(core_detail[Dataset %in% heat_ids], GeneSymbol ~ Dataset,
              value.var = "Statistic")
hm <- as.matrix(heat[, -1]); rownames(hm) <- heat$GeneSymbol
hm[hm > 3] <- 3; hm[hm < -3] <- -3
png(file.path(fig_dir, "Figure_core13_signed_statistic_heatmap.png"),
    width = 2800, height = 2200, res = 320)
pheatmap(hm, cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
         breaks = seq(-3, 3, length.out = 102), border_color = "grey85",
         fontsize = 9, main = "Thirteen shared OXPHOS leading-edge genes")
dev.off()

sens_plot <- score_tab[GeneSet == "GO:0006119" &
  (grepl("Primary|stratum|interaction|GAD|SVA", Analysis))]
sens_plot[, Label := paste(Cohort, Analysis, sep = "  ")]
p_sens <- ggplot(sens_plot, aes(Beta, reorder(Label, Beta), color = Cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.16) +
  geom_point(size = 2.4) +
  labs(x = "Adjusted difference in standardized OXPHOS score", y = NULL,
       title = "Within-cohort OXPHOS sensitivity analyses") +
  theme_minimal(base_size = 9) + theme(panel.grid.major.y = element_blank(),
                                   legend.position = "bottom")
ggsave(file.path(fig_dir, "Figure_OXPHOS_sensitivity_forest.png"), p_sens,
       width = 10, height = 8.5, dpi = 320)

# ------------------------------------------------------------------------------
# Standardized volcano atlases for all completed exploratory datasets
# ------------------------------------------------------------------------------
volcano_panel <- function(file, label, xlab = "Effect MDD minus control") {
  x <- fread(file)
  effect_col <- intersect(c("log2FoldChange", "logFC", "DeltaBeta"), names(x))[1]
  p_col <- intersect(c("pvalue", "P.Value"), names(x))[1]
  q_col <- intersect(c("padj", "adj.P.Val"), names(x))[1]
  gene_col <- intersect(c("GeneSymbol", "gene", "ProbeID"), names(x))[1]
  d <- data.table(effect = x[[effect_col]], p = x[[p_col]], q = x[[q_col]],
                  gene = if (!is.na(gene_col)) x[[gene_col]] else "")
  d <- d[is.finite(effect) & is.finite(p) & p > 0]
  d[, sig := !is.na(q) & q < 0.05]
  labs <- d[order(p)][seq_len(min(3, .N))]
  ggplot(d, aes(effect, -log10(p), color = sig)) +
    geom_point(size = 0.45, alpha = 0.55) +
    geom_text(data = labs, aes(label = gene), color = "black", size = 2.2,
              check_overlap = TRUE, vjust = -0.5, na.rm = TRUE) +
    scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey72")) +
    labs(title = label, subtitle = paste0("FDR < 0.05: ", sum(d$sig)),
         x = xlab, y = expression(-log[10](italic(P)))) +
    theme_minimal(base_size = 8) + theme(legend.position = "none",
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7))
}

blood_files <- c(
  GSE98793 = file.path(repo_root, "transcriptomics/blood/GSE98793/results_fixed/GSE98793_all_DEG_results.csv"),
  GSE39653 = file.path(repo_root, "transcriptomics/blood/GSE39653/results_fixed/GSE39653_all_DEG_results.csv"),
  GSE19738 = file.path(repo_root, "analysis_framework/GSE19738/results/deg_GSE19738.csv"),
  GSE251778 = file.path(repo_root, "transcriptomics/blood/GSE251778/results_fixed/GSE251778_DEG_overall.csv"))
brain_files <- c(
  GSE53987_Hippocampus = file.path(repo_root, "transcriptomics/brain/GSE53987/results_fixed/GSE53987_all_DE_results_hippocampus.csv"),
  GSE53987_BA46 = file.path(repo_root, "transcriptomics/brain/GSE53987/results_fixed/GSE53987_all_DE_results_Pre-frontal_cortex_BA46.csv"),
  GSE53987_Striatum = file.path(repo_root, "transcriptomics/brain/GSE53987/results_fixed/GSE53987_all_DE_results_Associative_striatum.csv"),
  GSE54563_ACC = file.path(repo_root, "transcriptomics/brain/GSE54563/results_fixed/GSE54563_All_DEG_Results.csv"),
  GSE80655_AnCg = file.path(repo_root, "transcriptomics/brain/GSE80655/results_fixed/GSE80655_all_DE_results_AnCg.csv"),
  GSE80655_DLPFC = file.path(repo_root, "transcriptomics/brain/GSE80655/results_fixed/GSE80655_all_DE_results_DLPFC.csv"),
  GSE80655_nAcc = file.path(repo_root, "transcriptomics/brain/GSE80655/results_fixed/GSE80655_all_DE_results_nAcc.csv"))
meth_files <- c(
  GSE88890 = file.path(repo_root, "transcriptomics/methylation/GSE88890/results_fixed/GSE88890_DMP_main_blockedSubject.csv"),
  GSE113725 = file.path(repo_root, "transcriptomics/methylation/GSE113725/results_fixed/GSE113725_DMP_all.csv"),
  GSE125105 = file.path(repo_root, "transcriptomics/methylation/GSE125105/results_fixed/GSE125105_DMP_all_covariates.csv"))

save_atlas <- function(files, outfile, ncol, xlab) {
  ps <- Map(function(f, nm) volcano_panel(f, nm, xlab), files, names(files))
  p <- patchwork::wrap_plots(ps, ncol = ncol)
  ggsave(file.path(fig_dir, outfile), p, width = 10.5,
         height = ceiling(length(ps) / ncol) * 3.5, dpi = 300)
}
save_atlas(blood_files, "Supplementary_volcano_blood.png", 2,
           "Expression effect MDD minus control")
save_atlas(brain_files, "Supplementary_volcano_brain.png", 2,
           "Expression effect MDD minus control")
save_atlas(meth_files, "Supplementary_volcano_methylation.png", 2,
           "Methylation effect MDD minus control")

sex_files <- c(
  GSE19738_female = file.path(repo_root, "analysis_framework/GSE19738/results/deg_GSE19738_female.csv"),
  GSE19738_male = file.path(repo_root, "analysis_framework/GSE19738/results/deg_GSE19738_male.csv"),
  GSE251778_female = file.path(repo_root, "transcriptomics/blood/GSE251778/results_fixed/GSE251778_DEG_female.csv"),
  GSE251778_male = file.path(repo_root, "transcriptomics/blood/GSE251778/results_fixed/GSE251778_DEG_male.csv"))
save_atlas(sex_files, "Supplementary_volcano_sex_stratified.png", 2,
           "Expression effect MDD minus control")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("Strict sensitivity and cross-cohort integration complete.\n")
cat("Outputs:", out_dir, "\n")

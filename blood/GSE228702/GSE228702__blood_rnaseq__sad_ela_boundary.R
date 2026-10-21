#!/usr/bin/env Rscript

# GSE228702: pre-defined OXPHOS analysis in a 2x2 peripheral-blood RNA-seq
# design (social anxiety disorder [SAD] x early-life adversity [ELA]).
# The deposited expression matrix was adjusted by the original investigators
# for estimated cell counts; age and gender are additionally included here.
#
# Run from the repository root:
#   Rscript analysis_framework/13_GSE228702_SAD_ELA_OXPHOS.R

repo_root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo_root, "analysis_framework", "_libpaths.R"))

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
  library(limma)
  library(fgsea)
  library(singscore)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(ggrepel)
  library(sandwich)
})

set.seed(228702)
base_dir <- file.path(repo_root, "analysis_framework", "GSE228702")
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "figures")
log_dir <- file.path(base_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(data_dir,
  "GSE228702_adjusted_expression_Cellcounts_granulatorAbis0_nnls.csv.gz")
series_file <- file.path(data_dir, "GSE228702_series_matrix.txt.gz")
if (!file.exists(matrix_file) || !file.exists(series_file)) {
  stop("Missing GSE228702 processed expression or series matrix")
}

cat("===== GSE228702 SAD x ELA OXPHOS validation =====\n")
gse <- getGEO(filename = series_file, getGPL = FALSE)
pd <- pData(gse)

strip_key <- function(x) trimws(sub("^[^:]+:", "", as.character(x)))
sample_id <- sub(",.*", "", pd$title)
meta <- data.frame(
  Sample = sample_id,
  GSM = pd$geo_accession,
  SAD = strip_key(pd$characteristics_ch1),
  ELA = strip_key(pd$characteristics_ch1.1),
  Age = suppressWarnings(as.numeric(strip_key(pd$characteristics_ch1.2))),
  Gender = strip_key(pd$characteristics_ch1.3),
  stringsAsFactors = FALSE)
meta$SAD <- factor(meta$SAD, levels = c("noSAD", "SAD"))
meta$ELA <- factor(meta$ELA, levels = c("noELA", "ELA"))
meta$Gender <- factor(tolower(meta$Gender), levels = c("female", "male"))
meta$Group <- factor(paste(meta$SAD, meta$ELA, sep = "_"),
                     levels = c("noSAD_noELA", "noSAD_ELA",
                                "SAD_noELA", "SAD_ELA"))
meta$Age_z <- as.numeric(scale(meta$Age))
rownames(meta) <- meta$Sample
cat("Metadata group counts:\n")
print(table(meta$Group, useNA = "ifany"))

x <- fread(cmd = paste("gzip -cd", shQuote(matrix_file)),
           showProgress = interactive())
gene_col <- names(x)[1]
genes <- as.character(x[[gene_col]])
x[, (gene_col) := NULL]
colnames(x) <- sub("^X(?=[0-9])", "", colnames(x), perl = TRUE)
missing_meta <- setdiff(colnames(x), rownames(meta))
if (length(missing_meta)) {
  stop("Expression columns absent from GEO metadata: ",
       paste(missing_meta, collapse = ", "))
}
meta <- meta[colnames(x), , drop = FALSE]
expr <- as.matrix(x)
storage.mode(expr) <- "double"
rownames(expr) <- genes
rm(x)

# Deterministic outcome-independent duplicate handling: retain the row with the
# greatest across-sample variance for each submitted gene symbol.
valid_gene <- !is.na(rownames(expr)) & nzchar(rownames(expr)) &
  rowMeans(is.finite(expr)) >= 0.95
expr <- expr[valid_gene, , drop = FALSE]
gene_var <- apply(expr, 1, var, na.rm = TRUE)
ord <- order(rownames(expr), -gene_var, na.last = TRUE)
expr <- expr[ord, , drop = FALSE]
expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
expr <- expr[apply(expr, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
cat("Expression:", nrow(expr), "unique genes x", ncol(expr), "samples\n")

# The deposited values are DESeq2-normalized counts after cell-composition
# adjustment and can be negative; treating them as log-expression creates
# extreme, heteroskedastic effect estimates.  For gene-level inference, apply a
# phenotype-blind rank-based inverse-normal transform within each gene.  This
# makes effects comparable and robust to a few very highly expressed rows while
# preserving each gene's between-person ordering.  Within-sample singscore below
# continues to use the original adjusted values.
rank_int <- function(z) {
  ok <- is.finite(z)
  ans <- rep(NA_real_, length(z))
  ans[ok] <- qnorm((rank(z[ok], ties.method = "average") - 0.5) / sum(ok))
  ans
}
expr_model <- t(apply(expr, 1, rank_int))
rownames(expr_model) <- rownames(expr)
colnames(expr_model) <- colnames(expr)

complete <- complete.cases(meta[, c("Group", "Age_z", "Gender")])
meta_cc <- droplevels(meta[complete, , drop = FALSE])
expr_score_cc <- expr[, rownames(meta_cc), drop = FALSE]
expr_cc <- expr_model[, rownames(meta_cc), drop = FALSE]
write.csv(meta_cc, file.path(data_dir, "GSE228702_metadata_parsed.csv"),
          row.names = FALSE)

# Four cell means make the marginal and interaction contrasts explicit.
design <- model.matrix(~ 0 + Group + Age_z + Gender, data = meta_cc)
colnames(design) <- sub("^Group", "", colnames(design))
if (qr(design)$rank != ncol(design)) stop("Rank-deficient GSE228702 design")
contrasts <- makeContrasts(
  SAD_marginal = (SAD_noELA + SAD_ELA - noSAD_noELA - noSAD_ELA) / 2,
  ELA_marginal = (noSAD_ELA + SAD_ELA - noSAD_noELA - SAD_noELA) / 2,
  SADxELA_interaction = SAD_ELA - SAD_noELA - noSAD_ELA + noSAD_noELA,
  levels = design)

fit <- lmFit(expr_cc, design)
fit_c <- eBayes(contrasts.fit(fit, contrasts), robust = TRUE, trend = TRUE)

go_tab <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = "GO:0006119", keytype = "GOALL",
  columns = c("SYMBOL", "ONTOLOGYALL")))
ox_genes <- sort(unique(go_tab$SYMBOL[go_tab$ONTOLOGYALL == "BP" &
                                       !is.na(go_tab$SYMBOL)]))
ox_present <- intersect(ox_genes, rownames(expr_cc))
ox_nuclear <- ox_present[!grepl("^MT-", ox_present)]
core13 <- intersect(c("CHCHD2", "COX4I1", "COX5A", "COX5B", "COX6C",
                      "NDUFA13", "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFC2",
                      "NDUFS3", "SIRT3", "UQCR11"), rownames(expr_cc))
writeLines(ox_present, file.path(out_dir, "GSE228702_GO0006119_genes_present.txt"))

# Sample-level rank scores.  These are independent of the phenotype labels.
ranked <- rankGenes(expr_score_cc)
meta_cc$Score_OXPHOS <- simpleScore(ranked, upSet = ox_present,
                                    centerScore = TRUE)$TotalScore
meta_cc$Score_nuclear <- simpleScore(ranked, upSet = ox_nuclear,
                                     centerScore = TRUE)$TotalScore
meta_cc$Score_core13 <- simpleScore(ranked, upSet = core13,
                                    centerScore = TRUE)$TotalScore
meta_cc$Score_OXPHOS_z <- as.numeric(scale(meta_cc$Score_OXPHOS))
meta_cc$Score_nuclear_z <- as.numeric(scale(meta_cc$Score_nuclear))
meta_cc$Score_core13_z <- as.numeric(scale(meta_cc$Score_core13))

contrast_estimate <- function(outcome, contrast_name, gene_set, caveat = "") {
  f <- as.formula(paste(outcome, "~ 0 + Group + Age_z + Gender"))
  m <- lm(f, data = meta_cc)
  V <- sandwich::vcovHC(m, type = "HC3")
  L <- rep(0, length(coef(m))); names(L) <- names(coef(m))
  group_names <- paste0("Group", c("noSAD_noELA", "noSAD_ELA",
                                    "SAD_noELA", "SAD_ELA"))
  if (contrast_name == "SAD_marginal") {
    L[group_names] <- c(-0.5, -0.5, 0.5, 0.5)
  } else if (contrast_name == "ELA_marginal") {
    L[group_names] <- c(-0.5, 0.5, -0.5, 0.5)
  } else if (contrast_name == "SADxELA_interaction") {
    L[group_names] <- c(1, -1, -1, 1)
  } else stop("Unknown contrast")
  b <- sum(L * coef(m)); se <- sqrt(drop(t(L) %*% V %*% L))
  z <- b / se
  data.table(Analysis = contrast_name, Outcome = gene_set,
             N = nobs(m), Beta = b, SE_HC3 = se,
             CI_low = b - 1.96 * se, CI_high = b + 1.96 * se,
             P = 2 * pt(abs(z), df = df.residual(m), lower.tail = FALSE),
             Caveat = caveat)
}

score_results <- rbindlist(lapply(c("SAD_marginal", "ELA_marginal",
                                    "SADxELA_interaction"), function(k) {
  rbind(
    contrast_estimate("Score_OXPHOS_z", k, "GO:0006119"),
    contrast_estimate("Score_nuclear_z", k, "Nuclear OXPHOS"),
    contrast_estimate("Score_core13_z", k, "Post-hoc core13",
                      "Discovery-derived score; descriptive only"))
}))
fwrite(score_results,
       file.path(out_dir, "GSE228702_OXPHOS_score_models.tsv"), sep = "\t")

# Gene-level tables plus two complementary rank-based pathway tests.
camera_results <- list()
fgsea_results <- list()
idx <- match(ox_present, rownames(expr_cc)); idx <- idx[!is.na(idx)]
for (nm in colnames(contrasts)) {
  de <- topTable(fit_c, coef = nm, number = Inf, sort.by = "P")
  de$GeneSymbol <- rownames(de)
  de <- de[, c("GeneSymbol", setdiff(colnames(de), "GeneSymbol"))]
  fwrite(as.data.table(de),
         file.path(out_dir, paste0("GSE228702_DE_", nm, ".tsv")), sep = "\t")

  cam <- camera(expr_cc, index = list(GO_0006119_OXPHOS = idx),
                design = design, contrast = contrasts[, nm])
  camera_results[[nm]] <- cbind(data.table(Analysis = nm),
                                as.data.table(cam, keep.rownames = "GeneSet"))

  stats <- fit_c$t[, nm]; names(stats) <- rownames(fit_c$t)
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  fg <- fgseaMultilevel(pathways = list(GO_0006119_OXPHOS = ox_present),
                       stats = stats, minSize = 10, maxSize = 500, eps = 0)
  fgsea_results[[nm]] <- cbind(data.table(Analysis = nm), as.data.table(fg))
}
fwrite(rbindlist(camera_results, fill = TRUE),
       file.path(out_dir, "GSE228702_OXPHOS_CAMERA.tsv"), sep = "\t")
fwrite(rbindlist(fgsea_results, fill = TRUE),
       file.path(out_dir, "GSE228702_OXPHOS_fgsea.tsv"), sep = "\t")

sample_out <- as.data.table(meta_cc, keep.rownames = "RowName")
fwrite(sample_out[, .(Sample, GSM, SAD, ELA, Group, Age, Gender,
                      Score_OXPHOS, Score_OXPHOS_z,
                      Score_nuclear_z, Score_core13_z)],
       file.path(out_dir, "GSE228702_sample_scores.tsv"), sep = "\t")

# Figures.
p_group <- ggplot(meta_cc, aes(Group, Score_OXPHOS_z, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.45, color = NA) +
  geom_boxplot(width = 0.20, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.08, alpha = 0.45, size = 1.2) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "GSE228702: OXPHOS score across SAD x ELA groups",
       subtitle = "Plot is unadjusted; inference uses age- and gender-adjusted HC3 contrasts",
       x = NULL, y = "Standardized OXPHOS singscore")
ggsave(file.path(fig_dir, "GSE228702_SAD_ELA_OXPHOS_groups.png"), p_group,
       width = 8.5, height = 6.4, dpi = 300)
ggsave(file.path(fig_dir, "GSE228702_SAD_ELA_OXPHOS_groups.pdf"), p_group,
       width = 8.5, height = 6.4)

forest <- score_results[Outcome == "GO:0006119"]
forest[, Analysis := factor(Analysis,
  levels = rev(c("SAD_marginal", "ELA_marginal", "SADxELA_interaction")))]
p_forest <- ggplot(forest, aes(Beta, Analysis)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.15,
                 color = "#24557A") +
  geom_point(size = 3, color = "#C0392B") +
  theme_bw(base_size = 12) +
  labs(title = "GSE228702 OXPHOS pre-specified contrasts",
       subtitle = "Adjusted standardized effects with HC3 95% confidence intervals",
       x = "Adjusted standardized effect", y = NULL)
ggsave(file.path(fig_dir, "GSE228702_OXPHOS_contrast_forest.png"), p_forest,
       width = 8.5, height = 5.6, dpi = 300)
ggsave(file.path(fig_dir, "GSE228702_OXPHOS_contrast_forest.pdf"), p_forest,
       width = 8.5, height = 5.6)

sad_de <- fread(file.path(out_dir, "GSE228702_DE_SAD_marginal.tsv"))
sad_de[, minus_log10_p := -log10(pmax(P.Value, .Machine$double.xmin))]
sad_de[, OXPHOS := GeneSymbol %in% ox_present]
sad_de[, Label := NA_character_]
top_ox_labels <- sad_de[OXPHOS == TRUE][order(P.Value)][1:min(12, .N), GeneSymbol]
sad_de[GeneSymbol %in% top_ox_labels, Label := GeneSymbol]
p_volcano <- ggplot(sad_de, aes(logFC, minus_log10_p, color = OXPHOS)) +
  geom_point(alpha = 0.42, size = 0.85) +
  scale_color_manual(values = c(`FALSE` = "grey65", `TRUE` = "#E67E22")) +
  geom_text_repel(aes(label = Label), na.rm = TRUE, max.overlaps = Inf,
                  min.segment.length = 0) +
  theme_bw(base_size = 12) +
  labs(title = "GSE228702 SAD marginal-effect transcriptome",
       subtitle = "Labelled genes are members of the pre-defined GO:0006119 set",
       x = "Adjusted expression contrast", y = expression(-log[10](P)))
ggsave(file.path(fig_dir, "GSE228702_SAD_marginal_volcano.png"), p_volcano,
       width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(fig_dir, "GSE228702_SAD_marginal_volcano.pdf"), p_volcano,
       width = 8.5, height = 6.5)

pc <- prcomp(t(expr_cc), scale. = FALSE)
pc_dat <- cbind(meta_cc, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
p_pca <- ggplot(pc_dat, aes(PC1, PC2, color = Group)) +
  geom_point(alpha = 0.75, size = 2) +
  theme_bw(base_size = 12) +
  labs(title = "GSE228702 adjusted-expression PCA",
       x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]))
ggsave(file.path(fig_dir, "GSE228702_expression_PCA.png"), p_pca,
       width = 8, height = 6.2, dpi = 300)
ggsave(file.path(fig_dir, "GSE228702_expression_PCA.pdf"), p_pca,
       width = 8, height = 6.2)

de_counts <- vapply(colnames(contrasts), function(nm) {
  d <- fread(file.path(out_dir, paste0("GSE228702_DE_", nm, ".tsv")))
  sum(d$adj.P.Val < 0.05, na.rm = TRUE)
}, integer(1))
summary_lines <- c(
  "GSE228702 SAD x ELA OXPHOS validation",
  "====================================",
  sprintf("Complete samples: %d", nrow(meta_cc)),
  paste("Group counts:", paste(names(table(meta_cc$Group)),
                               as.integer(table(meta_cc$Group)), collapse = "; ")),
  sprintf("Unique genes tested: %d", nrow(expr_cc)),
  sprintf("GO:0006119 genes present: %d", length(ox_present)),
  paste("FDR-significant genes:", paste(names(de_counts), de_counts, collapse = "; ")),
  "Primary interpretation must jointly consider singscore, CAMERA and fgsea.",
  "The deposited expression matrix was already adjusted for estimated cell counts.",
  "Gene-level inference used a phenotype-blind within-gene inverse-normal transform because adjusted counts include negative values.")
writeLines(summary_lines, file.path(out_dir, "GSE228702_analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

#!/usr/bin/env Rscript

# Strict processed-data sensitivity for GSE251780/GSE251778 using EpiDISH
# blood-cell estimates. Requires the outputs of script 15.

repo_root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo_root, "analysis_framework", "_libpaths.R"))
suppressPackageStartupMessages({
  library(data.table); library(limma); library(edgeR); library(EpiDISH)
  library(org.Hs.eg.db); library(missMethyl); library(sandwich); library(lmtest)
  library(ggplot2)
})
set.seed(25178016)

base_dir <- file.path(repo_root, "analysis_framework", "GSE251780")
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "figures")
matrix_file <- file.path(data_dir,
  "GSE251780_Mokhtari_MatrixProcessed_AverageAndPvalue.txt.gz")
metadata_file <- file.path(data_dir, "GSE251780_metadata.csv")
dmp_file <- file.path(out_dir, "GSE251780_MDD_DMP_all.tsv")
paired_file <- file.path(out_dir, "GSE251780_paired_sample_scores.tsv")
rna_count_file <- file.path(repo_root, "transcriptomics", "blood", "GSE251778",
  "data", "GSE251778_Mokhtari_rna_counts.csv.gz")
stopifnot(file.exists(matrix_file), file.exists(metadata_file), file.exists(dmp_file),
          file.exists(paired_file), file.exists(rna_count_file))

cat("===== GSE251780 EpiDISH cell-adjusted sensitivity =====\n")
keep_ids <- fread(dmp_file, select = "Probe")$Probe
gzcmd <- paste("gzip -cd", shQuote(matrix_file))
hdr <- names(fread(cmd = gzcmd, nrows = 0, check.names = FALSE))
beta_cols <- hdr[!grepl("Detection\\.Pval$", hdr) & hdr != hdr[1]]
beta_dt <- fread(cmd = gzcmd, select = c(hdr[1], beta_cols),
                 showProgress = interactive(), check.names = FALSE)
probe <- as.character(beta_dt[[1]])
ii <- match(keep_ids, probe)
if (anyNA(ii)) stop("QC-retained probes are missing")
beta <- as.matrix(beta_dt[ii, ..beta_cols]); storage.mode(beta) <- "double"
rownames(beta) <- keep_ids; colnames(beta) <- beta_cols
rm(beta_dt, probe, ii, keep_ids); gc()
beta[beta <= 0] <- 1e-6; beta[beta >= 1] <- 1 - 1e-6

data(centDHSbloodDMC.m, package = "EpiDISH")
ref_probe <- intersect(rownames(centDHSbloodDMC.m), rownames(beta))
if (length(ref_probe) < 250) stop("Too few EpiDISH reference probes")
cell <- EpiDISH::epidish(beta.m = beta[ref_probe, , drop = FALSE],
                         ref.m = centDHSbloodDMC.m[ref_probe, , drop = FALSE],
                         method = "RPC")$estF
cell <- as.data.frame(cell)
cell$SampleID <- rownames(cell)
fwrite(as.data.table(cell), file.path(out_dir, "GSE251780_EpiDISH_cell_fractions.tsv"),
       sep = "\t")

meta <- fread(metadata_file); meta[, SampleID := as.character(SampleID)]
meta <- meta[match(colnames(beta), SampleID)]
meta[, Group := factor(Group, levels = c("Control", "MDD"))]
meta[, Sex := factor(Sex, levels = c("Female", "Male"))]
cell <- cell[match(meta$SampleID, cell$SampleID), ]
cell_covars <- setdiff(colnames(cell), c("SampleID", "Neutro"))
model_dat <- cbind(as.data.frame(meta),
                   cell[, setdiff(colnames(cell), "SampleID"), drop = FALSE])

# Neutrophils are the omitted compositional reference cell type.
form <- as.formula(paste("~ Group + Sex +", paste(cell_covars, collapse = " + ")))
design <- model.matrix(form, data = model_dat)
M <- log2(beta / (1 - beta))
fit <- eBayes(lmFit(M, design), robust = TRUE, trend = TRUE)
dmp <- topTable(fit, coef = "GroupMDD", number = Inf, sort.by = "P")
dmp$Probe <- rownames(dmp)
fwrite(dmp, file.path(out_dir, "GSE251780_MDD_DMP_cell_adjusted_all.tsv"), sep = "\t")
fwrite(dmp[dmp$adj.P.Val < 0.05, ],
       file.path(out_dir, "GSE251780_MDD_DMP_cell_adjusted_FDR05.tsv"), sep = "\t")

ox_probes <- scan(file.path(out_dir, "GSE251780_OXPHOS_promoter_probes.txt"),
                  what = character(), quiet = TRUE)
ox_probes <- intersect(ox_probes, rownames(M))
ox_score <- as.numeric(scale(colMeans(t(scale(t(M[ox_probes, , drop = FALSE]))),
                                     na.rm = TRUE)))
model_dat$DNAm_OXPHOS_promoter_z <- ox_score
score_form <- as.formula(paste(
  "DNAm_OXPHOS_promoter_z ~ Group + Sex +", paste(cell_covars, collapse = " + ")))
score_fit <- lm(score_form, data = model_dat)
score_ct <- coeftest(score_fit, vcov. = vcovHC(score_fit, type = "HC3"))
b <- unname(score_ct["GroupMDD", 1]); se <- unname(score_ct["GroupMDD", 2])
fwrite(data.table(
  Analysis = "Cell-adjusted MDD effect on OXPHOS promoter DNAm score",
  N = nobs(score_fit), N_CpGs = length(ox_probes), Beta_MDD = b,
  SE_HC3 = se, CI_low = b - 1.96 * se, CI_high = b + 1.96 * se,
  P = unname(score_ct["GroupMDD", 4]),
  Cell_method = "EpiDISH RPC; neutrophil omitted as compositional reference"),
  file.path(out_dir, "GSE251780_OXPHOS_DNAm_score_cell_adjusted.tsv"), sep = "\t")

# Cell-adjusted paired pathway association.
paired <- fread(paired_file)
paired[, SampleID := as.character(SampleID)]
pd <- merge(paired[, .(SampleID, RNA_OXPHOS_all)],
            as.data.table(model_dat), by = "SampleID")
pair_form <- as.formula(paste(
  "RNA_OXPHOS_all ~ DNAm_OXPHOS_promoter_z + Group + Sex +",
  paste(cell_covars, collapse = " + ")))
pair_fit <- lm(pair_form, data = pd)
pair_ct <- coeftest(pair_fit, vcov. = vcovHC(pair_fit, type = "HC3"))
pb <- unname(pair_ct["DNAm_OXPHOS_promoter_z", 1])
pse <- unname(pair_ct["DNAm_OXPHOS_promoter_z", 2])
fwrite(data.table(
  Analysis = "Cell-adjusted paired RNA vs promoter DNAm OXPHOS scores",
  N = nobs(pair_fit), Beta_DNAm = pb, SE_HC3 = pse,
  CI_low = pb - 1.96 * pse, CI_high = pb + 1.96 * pse,
  P = unname(pair_ct["DNAm_OXPHOS_promoter_z", 4])),
  file.path(out_dir, "GSE251780_paired_pathway_cell_adjusted.tsv"), sep = "\t")

# Sex-stratified pathway slopes (also avoids comparing the two confounded
# RNA sequencing platforms across sex).
sex_res <- rbindlist(lapply(levels(pd$Sex), function(s) {
  z <- droplevels(pd[pd$Sex == s, ])
  sf <- as.formula(paste(
    "RNA_OXPHOS_all ~ DNAm_OXPHOS_promoter_z + Group +",
    paste(cell_covars, collapse = " + ")))
  ff <- lm(sf, data = z)
  cc <- coeftest(ff, vcov. = vcovHC(ff, type = "HC3"))
  bb <- unname(cc["DNAm_OXPHOS_promoter_z", 1])
  ss <- unname(cc["DNAm_OXPHOS_promoter_z", 2])
  data.table(Sex = s, N = nobs(ff), Beta_DNAm = bb, SE_HC3 = ss,
             CI_low = bb - 1.96 * ss, CI_high = bb + 1.96 * ss,
             P = unname(cc["DNAm_OXPHOS_promoter_z", 4]))
}))
fwrite(sex_res, file.path(out_dir,
  "GSE251780_paired_pathway_cell_adjusted_by_sex.tsv"), sep = "\t")

# Cell-adjusted targeted check of the sole gene-level FDR association from
# script 15, DNAJC15.  It is not treated as a disease-associated discovery.
dna_probes <- c("cg13561081", "cg11679069", "cg14729962", "cg12504148",
                "cg12201399", "cg00131557", "cg05035143", "cg15988970",
                "cg09677945")
dna_probes <- intersect(dna_probes, rownames(M))
dna_m <- colMeans(M[dna_probes, , drop = FALSE], na.rm = TRUE)
cnt <- fread(cmd = paste("gzip -cd", shQuote(rna_count_file)), check.names = FALSE)
ens <- sub("\\..*$", "", as.character(cnt$gene_ID))
dna_ens <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = "DNAJC15", keytype = "SYMBOL", columns = "ENSEMBL"))$ENSEMBL
jj <- which(ens %in% dna_ens)
ccols <- grep("_rna$", names(cnt), value = TRUE)
cmat <- as.matrix(cnt[jj, ..ccols]); storage.mode(cmat) <- "double"
rownames(cmat) <- ens[jj]; colnames(cmat) <- sub("_rna$", "", ccols)
all_counts <- as.matrix(cnt[, ..ccols]); storage.mode(all_counts) <- "double"
lib_size <- colSums(all_counts, na.rm = TRUE)
# Use total-library sizes. Calling cpm() on the one-gene subset would make the
# gene its own library denominator and incorrectly yield a constant vector.
dna_lcpm <- log2(sweep(cmat + 2, 2, lib_size + 4, "/") * 1e6)
dna_rna <- colMeans(dna_lcpm, na.rm = TRUE)
rm(all_counts, cnt, cmat, dna_lcpm); gc()
dna_dat <- as.data.table(model_dat)
dna_dat[, DNAJC15_DNAm := dna_m[match(SampleID, names(dna_m))]]
dna_dat[, DNAJC15_RNA := dna_rna[match(SampleID, names(dna_rna))]]
dna_form <- as.formula(paste(
  "DNAJC15_RNA ~ DNAJC15_DNAm + Group + Sex +", paste(cell_covars, collapse = " + ")))
dna_fit <- lm(dna_form, data = dna_dat)
dna_ct <- coeftest(dna_fit, vcov. = vcovHC(dna_fit, type = "HC3"))
db <- unname(dna_ct["DNAJC15_DNAm", 1]); dse <- unname(dna_ct["DNAJC15_DNAm", 2])
fwrite(data.table(Gene = "DNAJC15", N = nobs(dna_fit),
                  N_promoter_CpGs = length(dna_probes), Beta_DNAm = db,
                  SE_HC3 = dse, CI_low = db - 1.96 * dse,
                  CI_high = db + 1.96 * dse,
                  P = unname(dna_ct["DNAJC15_DNAm", 4])),
       file.path(out_dir, "GSE251780_DNAJC15_cell_adjusted.tsv"), sep = "\t")

# Cell-composition group differences, useful for diagnosing confounding.
cell_group <- rbindlist(lapply(c(colnames(centDHSbloodDMC.m)), function(ct) {
  dd <- data.frame(Value = model_dat[[ct]], Group = model_dat$Group,
                   Sex = model_dat$Sex)
  ff <- lm(Value ~ Group + Sex, data = dd)
  cct <- coeftest(ff, vcov. = vcovHC(ff, type = "HC3"))
  data.table(Cell = ct, Beta_MDD = unname(cct["GroupMDD", 1]),
             P = unname(cct["GroupMDD", 4]))
}))
cell_group[, FDR := p.adjust(P, "BH")]
fwrite(cell_group, file.path(out_dir, "GSE251780_EpiDISH_group_differences.tsv"),
       sep = "\t")

# Residualized visualization removes the large sex/platform offset.
rna_res <- residuals(lm(as.formula(paste(
  "RNA_OXPHOS_all ~ Group + Sex +", paste(cell_covars, collapse = " + "))), pd))
dna_res <- residuals(lm(as.formula(paste(
  "DNAm_OXPHOS_promoter_z ~ Group + Sex +", paste(cell_covars, collapse = " + "))), pd))
plot_dat <- data.frame(RNA_residual = rna_res, DNAm_residual = dna_res,
                       Group = pd$Group, Sex = pd$Sex)
p <- ggplot(plot_dat, aes(DNAm_residual, RNA_residual, color = Group, shape = Sex)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black") +
  scale_color_manual(values = c(Control = "#3B82B7", MDD = "#C74B50")) +
  theme_bw(base_size = 12) +
  labs(title = "Cell- and phenotype-adjusted paired OXPHOS scores",
       subtitle = "Residualized for diagnosis, sex and EpiDISH blood-cell fractions",
       x = "Residual promoter DNAm score", y = "Residual RNA score")
ggsave(file.path(fig_dir, "GSE251780_paired_DNAm_RNA_cell_adjusted.png"), p,
       width = 8, height = 6.2, dpi = 300)

summary_lines <- c(
  "GSE251780 EpiDISH strict sensitivity",
  "=====================================",
  sprintf("Samples: %d", nrow(model_dat)),
  sprintf("EpiDISH reference CpGs available: %d", length(ref_probe)),
  sprintf("Cell-adjusted FDR-significant MDD CpGs: %d",
          sum(dmp$adj.P.Val < 0.05, na.rm = TRUE)),
  "Model: MDD + sex + B + NK + CD4T + CD8T + Mono + Eosino",
  "Neutrophil is the omitted compositional reference.")
writeLines(summary_lines, file.path(out_dir,
  "GSE251780_cell_adjusted_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

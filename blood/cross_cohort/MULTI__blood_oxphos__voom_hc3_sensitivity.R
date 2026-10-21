# =============================================================================
# Step 3.1 corrections (2026-10-04)
# 1) GSE251778 CAMERA/FRY re-run with FULL voom object (weights), keep v$E
#    unweighted results side-by-side
# 2) GSE98793 nested models & GSE251780 paired models: OLS + hand-computed HC3
# 3) Sample overlap check: are GSE251780 paired samples a subset of GSE251778?
# 4) sessionInfo + corrected decision summary
# Does NOT overwrite step3 outputs.
# =============================================================================
suppressPackageStartupMessages({library(limma); library(edgeR)})

# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

OUT <- file.path(OUTPUT_ROOT, "cross_cohort/sensitivity_2026-10-04")
STRICT <- file.path(WORK_ROOT, "analysis_framework/results_strict")
D98793 <- file.path(WORK_ROOT, "transcriptomics/blood/GSE98793/results_fixed")
D251778 <- file.path(WORK_ROOT, "transcriptomics/blood/GSE251778")
D251780 <- file.path(WORK_ROOT, "analysis_framework/GSE251780")

log <- file(file.path(OUT, "step3.1_run.log"), open = "wt")
sink(log, split = TRUE)
on.exit(sink(), add = TRUE)
cat("=== Step 3.1 corrections ===\n"); cat(date(), "\n\n")

genes <- read.delim(file.path(OUT, "GO0006119_genes_frozen.tsv"),
                    stringsAsFactors = FALSE)
oxphos_all <- unique(genes$gene[genes$set == "all"])

# ---- HC3 helper -------------------------------------------------------------
hc3_table <- function(fit, term) {
  X <- model.matrix(fit); e <- residuals(fit); h <- hatvalues(fit)
  XtXi <- solve(crossprod(X))
  vcov_hc3 <- XtXi %*% t(X) %*% diag(e^2/(1-h)^2) %*% X %*% XtXi
  se_hc3 <- sqrt(diag(vcov_hc3))
  co <- summary(fit)$coefficients
  b <- co[term, "Estimate"]; df <- nrow(X) - ncol(X)
  se_ols <- co[term, "Std. Error"]; p_ols <- co[term, "Pr(>|t|)"]
  se3 <- se_hc3[term]; t3 <- b/se3
  data.frame(
    Beta = b,
    SE_OLS = se_ols, P_OLS = p_ols,
    CI_OLS_low = b - qt(0.975, df)*se_ols,
    CI_OLS_high = b + qt(0.975, df)*se_ols,
    SE_HC3 = se3, P_HC3 = 2*pt(abs(t3), df, lower.tail=FALSE),
    CI_HC3_low = b - qt(0.975, df)*se3,
    CI_HC3_high = b + qt(0.975, df)*se3,
    N = nrow(X), R2 = summary(fit)$r.squared,
    stringsAsFactors = FALSE)
}

# =============================================================================
# A. GSE98793 nested models: OLS + HC3
# =============================================================================
cat("################ A. GSE98793 OLS + HC3 ################\n")
sc <- read.delim(file.path(STRICT, "blood_OXPHOS_sample_scores.tsv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
d9 <- sc[sc$Cohort == "GSE98793", ]
d9$Group <- factor(d9$Group, levels = c("Control","MDD"))
d9$Sex <- factor(d9$Sex); d9$Batch <- factor(d9$Batch)

rows <- list()
for (sn in c("Score_all","Score_nuclear","Score_core13")) {
  d <- d9; d$y <- d[[sn]]
  for (mn in c("M0_unadj","M1_age_sex","M2_age_sex_batch")) {
    f <- switch(mn,
      M0_unadj = lm(y ~ Group, d),
      M1_age_sex = lm(y ~ Group + Age + Sex, d),
      M2_age_sex_batch = lm(y ~ Group + Age + Sex + Batch, d))
    r <- hc3_table(f, "GroupMDD")
    r$Score <- sn; r$Model <- mn
    rows[[length(rows)+1]] <- r
  }
}
tab9 <- do.call(rbind, rows)
tab9 <- tab9[, c("Score","Model","Beta","SE_OLS","CI_OLS_low","CI_OLS_high","P_OLS",
                 "SE_HC3","CI_HC3_low","CI_HC3_high","P_HC3","N","R2")]
write.table(tab9, file.path(OUT,"GSE98793_nested_models_OLS_HC3.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, na="")
print(tab9, digits=3)

# =============================================================================
# B. GSE251778: voom FULL object with weights (camera/fry)
# =============================================================================
cat("\n################ B. GSE251778 voom-weighted CAMERA/FRY ################\n")
cnt <- read.csv(gzfile(file.path(D251778,"data/GSE251778_Mokhtari_rna_counts.csv.gz")),
                row.names = 1, check.names = FALSE)
ensg <- cnt$gene_ID; cnt <- cnt[, -1, drop=FALSE]
colnames(cnt) <- sub("_rna$","",colnames(cnt))
emap <- read.delim(file.path(OUT,"ENSG_to_SYMBOL.tsv"), stringsAsFactors=FALSE)
sym <- emap$symbol[match(ensg, emap$ensembl_id)]
cnt <- cnt[!is.na(sym),,drop=FALSE]; sym2 <- sym[!is.na(sym)]
vars <- apply(cnt,1,var); o <- order(sym2,-vars); cnt <- cnt[o,]; sym2 <- sym2[o]
keep <- !duplicated(sym2); cnt <- cnt[keep,,drop=FALSE]; rownames(cnt) <- sym2[keep]

meta78 <- read.csv(file.path(D251778,"results_fixed/GSE251778_metadata_aligned.csv"),
                   stringsAsFactors=FALSE)
m <- meta78[match(colnames(cnt), as.character(meta78$SampleID)),]
m$Group <- factor(m$Group, levels=c("Control","MDD"))
m$Sex <- factor(m$Sex)
design <- model.matrix(~ Group + Sex, data=m)
y <- DGEList(counts=as.matrix(cnt))
keep2 <- filterByExpr(y, design=design); y <- y[keep2,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
v <- voom(y, design, plot=FALSE)
cat("voom genes:", nrow(v$E), " weights present:", !is.null(v$weights), "\n")

idx <- list(OXPHOS = which(rownames(v$E) %in% oxphos_all))
cat("OXPHOS genes in voom:", length(idx[[1]]), "\n")
cid <- match("GroupMDD", colnames(design))

# Weighted (full voom object)
cam_w_fixed <- camera(v, idx, design=design, contrast=cid, inter.gene.cor=0.01)
cam_w_est   <- camera(v, idx, design=design, contrast=cid, inter.gene.cor=NA)
fry_w       <- fry(v, idx, design=design, contrast=cid)
# Unweighted (v$E only) for side-by-side
cam_uw_fixed <- camera(v$E, idx, design=design, contrast=cid, inter.gene.cor=0.01)
cam_uw_est   <- camera(v$E, idx, design=design, contrast=cid, inter.gene.cor=NA)
fry_uw       <- fry(v$E, idx, design=design, contrast=cid)

cat("\n--- WEIGHTED (full voom object) ---\n")
print(cam_w_fixed); print(cam_w_est); print(fry_w)
cat("\n--- UNWEIGHTED (v$E only, kept for comparison) ---\n")
print(cam_uw_fixed); print(cam_uw_est); print(fry_uw)

mk <- function(x, label) {
  if (is.null(x)) return(data.frame(Test=label,NGenes=NA,Correlation=NA,
                                    Direction=NA,PValue=NA,stringsAsFactors=FALSE))
  data.frame(Test=label, NGenes=x$NGenes,
             Correlation=if(is.null(x$Correlation)) NA else x$Correlation,
             Direction=x$Direction, PValue=x$PValue, stringsAsFactors=FALSE)
}
b_out <- rbind(
  mk(cam_w_fixed, "WEIGHTED_CAMERA_fixed_r0.01"),
  mk(cam_w_est,   "WEIGHTED_CAMERA_estimated"),
  mk(fry_w,       "WEIGHTED_FRY"),
  mk(cam_uw_fixed,"UNWEIGHTED_CAMERA_fixed_r0.01"),
  mk(cam_uw_est,  "UNWEIGHTED_CAMERA_estimated"),
  mk(fry_uw,      "UNWEIGHTED_FRY"))
write.table(b_out, file.path(OUT,"GSE251778_CAMERA_FRY_voom_weighted.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, na="")

# =============================================================================
# C. Sample overlap: is GSE251780 paired subset of GSE251778?
# =============================================================================
cat("\n################ C. Sample overlap ################\n")
ps <- read.delim(file.path(D251780,"results/GSE251780_paired_sample_scores.tsv"),
                 stringsAsFactors=FALSE)
s78 <- as.character(meta78$SampleID)
s80 <- as.character(ps$SampleID)
ov <- intersect(s78, s80)
cat("GSE251778 SampleID n =", length(s78), "\n")
cat("GSE251780 paired SampleID n =", length(s80), "\n")
cat("overlap n =", length(ov), "\n")
cat("GSE251780 inside GSE251778? ", all(s80 %in% s78), "\n")
cat("overlap SampleIDs (first 10):", paste(head(ov,10), collapse=", "), "\n")

# =============================================================================
# D. GSE251780 paired RNA: OLS + HC3
# =============================================================================
cat("\n################ D. GSE251780 OLS + HC3 ################\n")
cf <- read.delim(file.path(D251780,"results/GSE251780_EpiDISH_cell_fractions.tsv"),
                stringsAsFactors=FALSE)
d <- merge(ps, cf, by="SampleID", all.x=TRUE, sort=FALSE)
d <- d[!is.na(d$B),]
d$Group <- factor(d$Group, levels=c("Control","MDD")); d$Sex <- factor(d$Sex)
rows2 <- list()
for (sn in c("RNA_OXPHOS_all","RNA_OXPHOS_core13")) {
  dd <- d; dd$y <- dd[[sn]]
  a <- lm(y ~ Group + Sex, dd)
  b <- lm(y ~ Group + Sex + B + NK + CD4T + CD8T + Mono + Eosino, dd)
  r1 <- hc3_table(a,"GroupMDD"); r1$Score<-sn; r1$Model<-"A_no_cell"
  r2 <- hc3_table(b,"GroupMDD"); r2$Score<-sn; r2$Model<-"B_cell_adjusted"
  rows2[[length(rows2)+1]] <- r1; rows2[[length(rows2)+1]] <- r2
}
tab80 <- do.call(rbind, rows2)
tab80 <- tab80[,c("Score","Model","Beta","SE_OLS","CI_OLS_low","CI_OLS_high","P_OLS",
                  "SE_HC3","CI_HC3_low","CI_HC3_high","P_HC3","N","R2")]
write.table(tab80, file.path(OUT,"GSE251780_paired_RNA_OLS_HC3.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, na="")
print(tab80, digits=3)

# ---- sessionInfo -----------------------------------------------------------
cat("\n################ sessionInfo ################\n")
print(sessionInfo())
cat("\n=== DONE ===\n"); cat(date(),"\n")
close(log)

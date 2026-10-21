#!/usr/bin/env Rscript
# Recalculate Plan A promoter score with v3.22.0 GO:0006119 EntrezIDs
# Preserves: beta matrix, sample labels, covariates, HC3 inference
# Changes: GO gene set definition (v3.21=151 -> v3.22=148)

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
})

# Step 1: Load EPIC annotation from existing prepared_mapping.rds (no need to reload EPIC package)
PROM_DIR <- "/Users/tan/Documents/Codex/2026-09-19/ple/outputs/secondary_planA_2026-10-05"
GSE_DIR  <- "/Users/tan/Documents/MDD-MultiOmics/analysis_framework/GSE251780"
OUT_DIR  <- file.path(PROM_DIR, "recalc_v322")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Step 1: Load existing mapping and annotation ===\n")
obj <- readRDS(file.path(PROM_DIR, "prepared_mapping.rds"))
ann <- obj$ann  # EPIC annotation with Probe, UCSC_RefGene_Name, UCSC_RefGene_Group
cat("EPIC annotation probes:", nrow(ann), "\n")

# Step 2: Use v3.22.0 org.Hs.eg.db (default .libPaths = system path = v3.22)
cat("\n=== Step 2: Query GO:0006119 with v3.22.0 ===\n")
# Force system path for org.Hs.eg.db
sys_libs <- .libPaths()
# Remove Developer path if present
dev_path <- "/Users/tan/Developer/envs/R/lib/R/library"
sys_libs_clean <- sys_libs[sys_libs != dev_path]
.libPaths(sys_libs_clean)

library(org.Hs.eg.db)
cat("org.Hs.eg.db version:", as.character(packageVersion("org.Hs.eg.db")), "\n")

# Query GO:0006119 EntrezIDs
go_v322 <- unique(na.omit(
  AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0006119", keytype = "GOALL", columns = "ENTREZID")$ENTREZID
))
cat("GO:0006119 EntrezIDs (v3.22):", length(go_v322), "\n")

# Step 3: Rebuild symbol-to-Entrez mapping (using v3.22)
cat("\n=== Step 3: Rebuild EPIC symbol->Entrez mapping (v3.22) ===\n")
syms <- unique(unlist(strsplit(ann$UCSC_RefGene_Name, ";", fixed = TRUE)))
syms <- syms[!is.na(syms) & nzchar(syms)]

canonical <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
known <- intersect(syms, canonical)
direct <- AnnotationDbi::select(org.Hs.eg.db, keys = known, keytype = "SYMBOL", columns = "ENTREZID")
names(direct)[1] <- "OldSymbol"
direct$MappingType <- "canonical"

unknown <- setdiff(syms, known)
aliases <- intersect(unknown, AnnotationDbi::keys(org.Hs.eg.db, keytype = "ALIAS"))
if (length(aliases) > 0) {
  alias <- AnnotationDbi::select(org.Hs.eg.db, keys = aliases, keytype = "ALIAS", columns = "ENTREZID")
  names(alias)[1] <- "OldSymbol"
  alias$MappingType <- "alias"
  mp <- rbind(direct, alias)
} else {
  mp <- direct
}
mp <- unique(mp[!is.na(mp$ENTREZID), ])
nk <- table(mp$OldSymbol)
mp$Ambiguous <- nk[mp$OldSymbol] > 1

safe <- mp[!mp$Ambiguous, ]
idmap <- setNames(safe$ENTREZID, safe$OldSymbol)
cat("Mapped symbols:", length(idmap), "\n")

# Step 4: Rebuild promoter mapping for GO:0006119 set
cat("\n=== Step 4: Rebuild promoter mapping (v3.22 GO targets) ===\n")
targets <- go_v322  # Entrez IDs from v3.22

rows <- lapply(seq_len(nrow(ann)), function(i) {
  sy <- strsplit(ann$UCSC_RefGene_Name[i], ";", fixed = TRUE)[[1]]
  gr <- strsplit(ann$UCSC_RefGene_Group[i], ";", fixed = TRUE)[[1]]
  if (length(sy) != length(gr)) return(NULL)
  ids <- unname(idmap[sy])
  ok <- gr %in% c("TSS1500", "TSS200", "1stExon") & !is.na(ids) & ids %in% targets
  if (!any(ok)) return(NULL)
  unique(data.frame(
    Probe = ann$Probe[i],
    EntrezID = ids[ok],
    OldSymbol = sy[ok],
    Region = gr[ok]
  ))
})
pm_v322 <- unique(do.call(rbind, rows))
cat("Promoter probes (v3.22 GO):", length(unique(pm_v322$Probe)), "\n")
cat("Promoter genes (v3.22 GO):", length(unique(pm_v322$EntrezID)), "\n")

fwrite(pm_v322, file.path(OUT_DIR, "promoter_mapping_v322.tsv"), sep = "\t")

# Step 5: Load beta matrix and metadata
cat("\n=== Step 5: Load beta matrix and covariates ===\n")
x <- fread(file.path(PROM_DIR, "target_promoter_beta.tsv"), check.names = FALSE)
meta <- fread(file.path(GSE_DIR, "data/GSE251780_metadata.csv"))
cf <- fread(file.path(GSE_DIR, "results/GSE251780_EpiDISH_cell_fractions.tsv"))

d <- merge(meta, cf, by = "SampleID", sort = FALSE)
d$SampleID <- as.character(d$SampleID)
d$Group <- factor(d$Group, levels = c("Control", "MDD"))
d$Sex <- factor(d$Sex)
stopifnot(nrow(d) == 147, !anyDuplicated(d$SampleID), all(d$SampleID %in% names(x)))

beta <- as.matrix(x[, d$SampleID, with = FALSE])
rownames(beta) <- x[[1]]
storage.mode(beta) <- "double"
stopifnot(all(is.finite(beta)), all(beta >= 0 & beta <= 1))

M <- log2(pmin(pmax(beta, 1e-6), 1 - 1e-6) / (1 - pmin(pmax(beta, 1e-6), 1 - 1e-6)))

# Step 6: Filter promoter mapping to probes present in beta matrix
cat("\n=== Step 6: Compute gene-level score (v3.22) ===\n")
cat("Probes in beta matrix:", nrow(beta), "\n")
cat("Probes in promoter mapping:", length(unique(pm_v322$Probe)), "\n")
pm_v322 <- pm_v322[pm_v322$Probe %in% rownames(M), ]
cat("Probes after filtering to beta matrix:", length(unique(pm_v322$Probe)), "\n")
cat("Genes after filtering:", length(unique(pm_v322$EntrezID)), "\n")
ids <- unique(pm_v322$EntrezID)
geneM <- t(vapply(ids, function(id) {
  colMeans(M[unique(pm_v322$Probe[pm_v322$EntrezID == id]), , drop = FALSE])
}, numeric(nrow(d))))
rownames(geneM) <- ids

sdg <- apply(geneM, 1, sd)
ok <- is.finite(sdg) & sdg > 0
geneM <- geneM[ok, , drop = FALSE]
cat("Genes with non-zero SD:", nrow(geneM), "\n")

Z <- t(scale(t(geneM)))
present <- rownames(Z)
score <- colMeans(Z[present, , drop = FALSE])
d$Score <- as.numeric(scale(score))

cat("Score computed. N =", nrow(d), "\n")
cat("MDD mean score:", round(mean(d$Score[d$Group == "MDD"]), 4), "\n")
cat("Control mean score:", round(mean(d$Score[d$Group == "Control"]), 4), "\n")

# Step 7: Fit models with diagnostics
cat("\n=== Step 7: Fit OLS/HC3 models ===\n")

run_model <- function(data, formula_str, label) {
  f <- as.formula(formula_str)
  fit <- lm(f, data)
  b <- coef(fit)["GroupMDD"]
  se <- tryCatch(
    sqrt(vcovHC(fit, type = "HC3")["GroupMDD", "GroupMDD"]),
    error = function(e) NA_real_
  )
  if (is.finite(se) && se > 0) {
    ci <- b + c(-1, 1) * qt(.975, df.residual(fit)) * se
    p_hc3 <- 2 * pt(-abs(b / se), df.residual(fit))
  } else {
    ci <- c(NA_real_, NA_real_)
    p_hc3 <- NA_real_
  }
  p_ols <- coef(summary(fit))["GroupMDD", "Pr(>|t|)"]
  data.frame(
    Model = label,
    N = nobs(fit),
    Beta = b,
    SE_HC3 = se,
    CI_low = ci[1],
    CI_high = ci[2],
    P_HC3 = p_hc3,
    P_OLS = p_ols,
    Max_leverage = max(hatvalues(fit)),
    stringsAsFactors = FALSE
  )
}

# Find sample 192 (max leverage from original analysis)
fit_full <- lm(Score ~ Group + Sex + B + NK + CD4T + CD8T + Mono + Eosino, d)
leverage_full <- hatvalues(fit_full)
max_lev_idx <- which.max(leverage_full)
sample_192 <- d$SampleID[max_lev_idx]
cat("Max leverage sample:", sample_192, " leverage =", round(leverage_full[max_lev_idx], 4), "\n")

results <- list()

# Model 1: Sex-adjusted
results[[1]] <- run_model(d, "Score ~ Group + Sex", "sex_adjusted")

# Model 2: Full cell-adjusted (original)
results[[2]] <- run_model(d, "Score ~ Group + Sex + B + NK + CD4T + CD8T + Mono + Eosino", "cell_adjusted_full")

# Model 3: Exclude sample 192, same covariates (HC3 may be non-estimable)
d_no192 <- d[d$SampleID != sample_192, ]
results[[3]] <- run_model(d_no192, "Score ~ Group + Sex + B + NK + CD4T + CD8T + Mono + Eosino", "cell_adjusted_drop192")

# Model 4: Full sample, remove Eosino
results[[4]] <- run_model(d, "Score ~ Group + Sex + B + NK + CD4T + CD8T + Mono", "cell_adjusted_noEosino")

# Model 5: Drop sample 192 AND remove Eosino
results[[5]] <- run_model(d_no192, "Score ~ Group + Sex + B + NK + CD4T + CD8T + Mono", "cell_adjusted_drop192_noEosino")

res_df <- do.call(rbind, results)
res_df$Set <- "GO0006119_v322"
res_df$Excluded_Sample <- ifelse(grepl("drop192", res_df$Model), sample_192, "none")

print(res_df)

# Step 8: Coverage summary
cat("\n=== Step 8: Coverage summary ===\n")
cov_df <- data.frame(
  Set = "GO0006119_v322",
  Defined_genes = length(go_v322),
  Covered_genes = length(present),
  Probes = length(unique(pm_v322$Probe)),
  stringsAsFactors = FALSE
)
print(cov_df)

# Step 9: Save outputs
fwrite(res_df, file.path(OUT_DIR, "PlanA_recalc_v322_results.tsv"), sep = "\t")
fwrite(cov_df, file.path(OUT_DIR, "PlanA_recalc_v322_coverage.tsv"), sep = "\t")
fwrite(d[, c("SampleID", "Score")], file.path(OUT_DIR, "PlanA_recalc_v322_sample_scores.tsv"), sep = "\t")

capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo_recalc_v322.txt"))

cat("\n=== DONE ===\n")
cat("Outputs saved to:", OUT_DIR, "\n")

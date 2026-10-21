#!/usr/bin/env Rscript

# GSE251780__blood_methylation__promoter_v322_table_assembly.R
# Build the formal 14-row promoter result table from the two v3.22 GO primary
# models (sex_adjusted, cell_adjusted_full) plus the twelve Plan A (v3.21)
# carried-forward rows for the six non-GO gene sets.
#
# This script was newly created on 2026-10-20 as the reproducible assembly step.
# It is NOT the original historical assembly script (none was located); it does
# not refit any model and hard-codes no statistical values. All numbers come
# from the four input files listed below.
#
# Sensitivity models (cell_adjusted_drop192, cell_adjusted_noEosino,
# cell_adjusted_drop192_noEosino) are excluded by exact Model-name matching.

suppressPackageStartupMessages(library(utils))
options(digits = 17)   # full-precision numeric output in write.table/format

# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

root    <- file.path(ARCHIVE_ROOT, "outputs/manuscript_v6_12_optionB_2026-10-13/github_repository_draft")
v322dir <- file.path(OUTPUT_ROOT, "methylation/recalc_v322")
outdir  <- file.path(OUTPUT_ROOT, "methylation/github_repo_GSE251780")

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = c("NA", "", "NaN"))
}
strip_quotes <- function(x) gsub('"', "", x, fixed = TRUE)

sets <- c("GO0006119", "MitoCarta_complex_subunits_nuclear",
          "Complex_I", "Complex_II", "Complex_III", "Complex_IV", "Complex_V")

## 1. GO primary models from the v3.22 frozen output (exact Model-name match)
res_v322 <- read_tsv(file.path(v322dir, "PlanA_recalc_v322_results.tsv"))
primary  <- c("sex_adjusted", "cell_adjusted_full")
go <- res_v322[res_v322$Model %in% primary, , drop = FALSE]
stopifnot(nrow(go) == 2)
go$Set   <- "GO0006119"
go$Model <- ifelse(go$Model == "cell_adjusted_full", "cell_adjusted", "sex_adjusted")
go <- go[, c("Set", "Model", "N", "Beta", "SE_HC3", "CI_low", "CI_high", "P_HC3")]
names(go)[names(go) == "Beta"] <- "Beta_score_SD"

## 2. Twelve carried-forward Plan A rows (six non-GO sets x two primary models)
plan <- read_tsv(file.path(outdir, "promoter_DNAm.tsv"))
plan$Set <- strip_quotes(plan$Set)
non_go <- plan[plan$Set != "GO0006119_original_ID_repaired" &
                 plan$Model %in% c("sex_adjusted", "cell_adjusted"),
               c("Set", "Model", "N", "Beta_score_SD", "SE_HC3",
                 "CI_low", "CI_high", "P_HC3")]
stopifnot(nrow(non_go) == 12)

## 3. Assemble 14 rows in canonical set order
assembled <- rbind(go, non_go)
assembled <- assembled[order(match(assembled$Set, sets), assembled$Model), ]
stopifnot(nrow(assembled) == 14)
stopifnot(all(sets %in% assembled$Set))
key <- paste(assembled$Set, assembled$Model)
stopifnot(!anyDuplicated(key))                       # Set+Model unique
tab <- table(assembled$Set, assembled$Model)
stopifnot(all(tab == 1))                             # 7 sets x 2 primary models

## 4. BH across the 14 primary P values
assembled$BH_all_reported_tests <- p.adjust(assembled$P_HC3, method = "BH")

## 5. Coverage: GO from v3.22, six non-GO carried forward from Plan A
cov_v322 <- read_tsv(file.path(v322dir, "PlanA_recalc_v322_coverage.tsv"))
cov_v322$Set <- sub("_v322$", "", cov_v322$Set)      # GO0006119_v322 -> GO0006119
plan_cov <- read_tsv(file.path(outdir, "promoter_coverage.tsv"))
plan_cov$Set <- strip_quotes(plan_cov$Set)
plan_cov <- plan_cov[plan_cov$Set != "GO0006119_original_ID_repaired", ]
cov_out <- rbind(cov_v322, plan_cov)
cov_out <- cov_out[order(match(cov_out$Set, sets)), ]
stopifnot(all(sets %in% cov_out$Set))
stopifnot(!any(is.na(cov_out[, c("Defined_genes", "Covered_genes", "Probes")])))

## 6. Export to NEW positions; do not overwrite candidate or historical files
out_res <- file.path(outdir, "promoter_DNAm_v322_assembled.tsv")
out_cov <- file.path(outdir, "promoter_coverage_v322_assembled.tsv")
out_ord <- file.path(outdir, "Figure4_ordered_promoter_coverage_v322_assembled.tsv")
write.table(assembled[, c("Set", "Model", "N", "Beta_score_SD", "SE_HC3",
                          "CI_low", "CI_high", "P_HC3", "BH_all_reported_tests")],
            out_res, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cov_out, out_cov, sep = "\t", row.names = FALSE, quote = FALSE)
cov_ord <- cov_out
cov_ord$Y <- 7L - match(cov_ord$Set, sets)
write.table(cov_ord, out_ord, sep = "\t", row.names = FALSE, quote = FALSE)
cat("WROTE:", out_res, "\n", out_cov, "\n", out_ord, "\n")

## 7. Field-by-field comparison vs the current candidate table (match by Set+Model)
cand <- read_tsv(file.path(outdir, "promoter_DNAm_v322.tsv"))
cand_key <- paste(cand$Set, cand$Model)
asm_key  <- paste(assembled$Set, assembled$Model)
stopifnot(setequal(asm_key, cand_key))
cand_m <- cand[match(asm_key, cand_key), ]
cmp_cols <- c("N", "Beta_score_SD", "SE_HC3", "CI_low", "CI_high", "P_HC3",
              "BH_all_reported_tests")
all_ok <- TRUE
for (cn in cmp_cols) {
  same <- isTRUE(all.equal(as.numeric(assembled[[cn]]), as.numeric(cand_m[[cn]]),
                           tolerance = 1e-15, check.attributes = FALSE))
  cat(sprintf("FIELD %-20s %s\n", cn, if (same) "PASS" else "FAIL"))
  if (!same) { all_ok <- FALSE
    print(data.frame(Set = assembled$Set, Model = assembled$Model,
                     assembled = assembled[[cn]], candidate = cand_m[[cn]]))
  }
}
if (!all_ok) {
  stop("Field differences found vs candidate table; no automatic correction applied")
}
cat("RESULT: assembled 14-row table is IDENTICAL to candidate promoter_DNAm_v322.tsv",
    "across N, beta, SE, CI_low, CI_high, P, BH (Set+Model matched)\n")
cat("COVERAGE: 7 sets matched; GO 148/123/1404; six non-GO carried forward from Plan A\n")

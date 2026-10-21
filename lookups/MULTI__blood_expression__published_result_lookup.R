# ==============================================================================
# 07_validate_blood_results.R
# PURPOSE: Reproducible cross-cohort validation of the two completed blood
#          transcriptome analyses.  It does not call a gene/pathway "replicated"
#          unless direction agrees and the validation cohort passes BH FDR < .05.
# INPUT  : fixed all-gene tables and fixed GO-BP GSEA tables for GSE98793 and
#          GSE39653.
# OUTPUT : gene-level and pathway-level replication tables plus a text summary.
# NOTE   : GO terms are highly redundant.  Counts of replicated GO terms must
#          be interpreted as biological themes, not independent discoveries.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    <- "<DATA_ROOT>"      # raw GEO / supplementary files (see README "Inputs")
WORK_ROOT    <- "<WORK_ROOT>"       # author analysis workspace (intermediate result tables)
ARCHIVE_ROOT <- "<ARCHIVE_ROOT>"    # author manuscript output archive
OUTPUT_ROOT  <- "./output"          # results are written here (relative to repository root)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_arg)) dirname(normalizePath(script_arg)) else getwd()
repo_root <- WORK_ROOT
out_dir <- file.path(OUTPUT_ROOT, "lookups")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

deg_98793_file <- file.path(repo_root, "transcriptomics/blood/GSE98793/results_fixed/GSE98793_all_DEG_results.csv")
deg_39653_file <- file.path(repo_root, "transcriptomics/blood/GSE39653/results_fixed/GSE39653_all_DEG_results.csv")
gsea_98793_file <- file.path(repo_root, "transcriptomics/blood/GSE98793/results_fixed/GSE98793_GSEA_BP.csv")
gsea_39653_file <- file.path(repo_root, "transcriptomics/blood/GSE39653/results_fixed/GSE39653_GSEA_BP.csv")

required <- c(deg_98793_file, deg_39653_file, gsea_98793_file, gsea_39653_file)
if (any(!file.exists(required))) {
  stop("Missing required fixed result(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

# Gene-level validation: discovery set is all GSE98793 genes at BH FDR < .05.
a <- fread(deg_98793_file)
b <- fread(deg_39653_file)
a <- a[adj.P.Val < 0.05, .(GeneSymbol, logFC_98793 = logFC,
                            P_98793 = P.Value, FDR_98793 = adj.P.Val)]
b <- b[, .(GeneSymbol, logFC_39653 = logFC,
           P_39653 = P.Value, FDR_39653 = adj.P.Val)]
gene_rep <- merge(a, b, by = "GeneSymbol", all.x = TRUE)
gene_rep[, same_direction := sign(logFC_98793) == sign(logFC_39653)]
gene_rep[, replication_FDR_within_discovery_set :=
           p.adjust(ifelse(same_direction, P_39653, 1), method = "BH")]
gene_rep[, replicated := same_direction & replication_FDR_within_discovery_set < 0.05]
setorder(gene_rep, replication_FDR_within_discovery_set)
fwrite(gene_rep, file.path(out_dir, "blood_gene_replication_GSE98793_in_GSE39653.tsv"), sep = "\t")

# Pathway validation: require same NES direction and BH FDR < .05 in both.
ga <- fread(gsea_98793_file)
gb <- fread(gsea_39653_file)
keep <- c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust")
ga <- ga[, ..keep]
gb <- gb[, ..keep]
setnames(ga, setdiff(keep, "ID"), paste0(setdiff(keep, "ID"), "_98793"))
setnames(gb, setdiff(keep, "ID"), paste0(setdiff(keep, "ID"), "_39653"))
path_rep <- merge(ga, gb, by = "ID")
path_rep[, discovery_sig := p.adjust_98793 < 0.05]
path_rep[, same_direction := sign(NES_98793) == sign(NES_39653)]
path_rep[, replicated := discovery_sig & same_direction & p.adjust_39653 < 0.05]
setorder(path_rep, p.adjust_39653)
fwrite(path_rep, file.path(out_dir, "blood_GO_BP_replication_GSE98793_in_GSE39653.tsv"), sep = "\t")

disc_path <- path_rep[discovery_sig == TRUE]
summary_lines <- c(
  "Blood cross-cohort validation summary",
  "=====================================",
  sprintf("GSE98793 discovery genes at BH FDR < .05: %d", nrow(a)),
  sprintf("Mapped in GSE39653: %d", sum(!is.na(gene_rep$P_39653))),
  sprintf("Same gene-level direction: %d", sum(gene_rep$same_direction, na.rm = TRUE)),
  sprintf("Genes replicated at validation-set BH FDR < .05: %d", sum(gene_rep$replicated, na.rm = TRUE)),
  sprintf("GSE98793 GO-BP terms at BH FDR < .05 and matched: %d", nrow(disc_path)),
  sprintf("Same pathway direction: %d (%.1f%%)", sum(disc_path$same_direction),
          100 * mean(disc_path$same_direction)),
  sprintf("GO-BP terms at BH FDR < .05 in both, same direction: %d",
          sum(disc_path$replicated)),
  sprintf("NES correlation across discovery-significant matched terms: %.3f",
          cor(disc_path$NES_98793, disc_path$NES_39653)),
  "Interpret GO counts as redundant themes, not independent findings."
)
writeLines(summary_lines, file.path(out_dir, "blood_cross_cohort_validation_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

# GSE39653 — Peripheral Blood Microarray (Discovery Screen, Cohort 2)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE39653](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39653) |
| **Title** (approx.) | Whole-blood gene expression in mood disorders vs controls (Savitz et al.) |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina HumanHT-12 V4.0 expression beadchip (GPL10558) |
| **Technology** | Microarray |
| **Reference** | Savitz J, Frank MB, Victor T, et al. *Brain Behav Immun*. 2013;31:161–171. doi:10.1016/j.bbi.2012.10.007 |

---

## 2. Sample Design

| Metric | Value |
|--------|-------|
| **Total analytic samples** | 45 (MDD and control only; bipolar disorder samples excluded per v6.12 scope) |
| **Genes (post-QC, retained)** | 19,948 |
| **Tissue** | Whole blood / PBMC |
| **Covariates available** | Diagnosis only in the central OXPHOS model |

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Second cohort contributing to the initial broad GO-BP screen; one of three microarray blood cohorts in the participant-score synthesis |
| **Manuscript sections** | §2.1–§2.4; Tables 1–3; Tables S2, S15–S17 |
| **Screening contribution** | Broad GO-BP returned 6,174 terms; 124 at BH q<0.05. OXPHOS NES = **−2.01, q = 0.000776**, concordant in direction with GSE98793 |

---

## 4. Processing in This Re-analysis

From Table S15:

- MDD/control samples only; positive expression in ≥10% samples; conditional log₂(x+1) where needed.
- Probe → symbol via `illuminaHumanv4.db`; duplicate symbols collapsed to highest-mean probe.
- Gene-level design: `diagnosis` only (no age/sex/batch available in the retained model); ordinary eBayes.
- OXPHOS pathway: `camera` fixed r = 0.01 and estimated r; directional FRY.
- Participant score: centered directional singscore, standardised within cohort; OLS with HC3; MCP-counter proxy adjustment.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE39653__blood_microarray__differential_expression.R` |
| Required expression/phenotype input | The script accepts the first available of `GSE39653_blood_data_manual.RData`, `GSE39653_clean_data.rds`, or `GSE39653_series_matrix.txt`; the current script lists workstation-specific candidate paths. These files are not included here. |
| Additional input | `reference/MCPcounter/genes.txt` is read if the MCP-counter marker analysis is available. It is not stored in the GSE-specific directory. |
| Main outputs written by the script | `GSE39653_expr_clean.rds`, `GSE39653_gene_expression.rds`, `GSE39653_metadata.rds`, `GSE39653_all_DEG_results.csv`, `GSE39653_sig_DEGs_adjP0.05.csv`, `GSE39653_immune_fullBH_stats.csv`, `GSE39653_GSEA_BP.csv`, volcano TIFF/PDF, and a report/session record. |
| Selected outputs included here | results/by_cohort/GSE39653/ *(results are not included in this scripts-only release; see manuscript archive)* contains the all/significant DEG tables, immune and GSEA tables, volcano PDF, and report. Intermediate RDS files are not included. |
| Path/run note | The script writes to an absolute project path under `transcriptomics/blood/GSE39653/results_fixed`; it is not portable without a path configuration change. |

---

## 5. Key OXPHOS Results in v6.12

| Test | Value |
|------|-------|
| Preranked GO-BP OXPHOS NES (discovery) | −2.01, q = 0.000776 (downward) |
| CAMERA fixed r = 0.01 (119 genes) | P = 8.15 × 10⁻⁴, within-method q = 0.00122 |
| CAMERA estimated correlation (r = 0.063) | P = 0.0955, q = 0.143 |
| Directional FRY | P = 0.0869, q = 0.130 |
| Rank-score β (diagnosis only) | Included in the three-cohort REML synthesis |
| Mean-z score base → +cell proxies | −0.145 → −0.187; adjusted nominal P = 0.0298, q = 0.195 across 8 tests |

**Interpretation**: As in GSE98793, the fixed-correlation pathway result did not survive estimated-correlation correction. The cell-marker-adjusted mean-z estimate reached nominal significance but not the eight-test family correction.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39653>
- Obtain deposited series matrix directly from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `AnnotationDbi, clusterProfiler, dplyr, ggplot2, illuminaHumanv4.db, limma, org.Hs.eg.db, tidyr`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE39653__blood_microarray__differential_expression.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Savitz J, Frank MB, Victor T, et al. Inflammation and neurological disease-related genes are differentially expressed in depressed patients with mood disorders and correlate with morphometric and functional imaging abnormalities. *Brain Behav Immun*. 2013;31:161–171.
- Wu D, Smyth GK. Camera. *Nucleic Acids Res*. 2012;40:e133.

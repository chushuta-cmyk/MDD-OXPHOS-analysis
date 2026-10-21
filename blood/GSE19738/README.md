# GSE19738 — Stimulated Whole-Blood Microarray (Basal Subset)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE19738](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19738) |
| **Title** | Stimulated gene expression profiles as a blood marker of major depressive disorder (non-medicated MDD vs controls) |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina HumanRef-8 v3.0 expression beadchip (GPL6848) |
| **Technology** | Microarray (whole-blood, ex-vivo stimulated with two conditions) |
| **Reference** | Spijker S, Van Zanten JS, De Jong S, et al. *Biol Psychiatry*. 2010;68(2):179–186. doi:10.1016/j.biopsych.2010.03.017 |

---

## 2. Sample Design

| Metric | Value |
|--------|-------|
| **GEO total samples** | 132 (basal + stimulated, MDD and control) |
| **Analytic samples in v6.12** | 67 (**basal condition only**; stimulated duplicates excluded) |
| **Genes (post-QC, retained)** | 9,510 (finite probes after filtering) |
| **Tissue** | Whole blood |
| **Covariates available** | Age, sex |

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Third microarray blood cohort in the OXPHOS score synthesis; basal expression only (per Spijker protocol) |
| **Manuscript sections** | §2.1–§2.4; Tables 1–3; Tables S2, S15–S17 |
| **Notable feature** | Sex and platform are **not** confounded here; this cohort contributes the strongest fixed-correlation CAMERA signal among the three microarrays |

---

## 4. Processing in This Re-analysis

From Table S15:

- Basal samples only; finite probes retained; GPL6848 symbols used; duplicate symbols collapsed to highest-mean probe.
- No stimulated duplicate samples carried forward.
- Gene-level design: `diagnosis + age + sex`; robust/trend eBayes.
- OXPHOS pathway: `camera` fixed r = 0.01 and estimated r; directional FRY.
- Participant score: centered directional singscore, standardised within cohort; OLS with HC3; MCP-counter proxy adjustment.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE19738__blood_microarray__differential_expression.R` |
| Required input | The script retrieves the GEO series matrix through GEOquery and, by default, caches `GSE19738_series_matrix.txt.gz` under `scripts/expression/GSE19738/cache/` (or an optional cache directory supplied on the command line). GPL annotation is also retrieved through GEOquery when needed. These inputs are not included in this repository. |
| Main outputs written by the script | Under `scripts/expression/GSE19738/results/`: `GSE19738_gene_expression.rds`, `GSE19738_metadata.rds`, `deg_GSE19738.csv`, `deg_GSE19738_female.csv`, `deg_GSE19738_male.csv`, `GSE19738_GSEA_GO_BP*.tsv`, `GSE19738_predefined_GO0006119*.tsv`, `GSE19738_analysis_summary.tsv`, volcano PDFs, and report/session files. Unified DEG TSVs are written under `scripts/expression/GSE19738/unified_tables/`. |
| Selected outputs included here | results/by_cohort/GSE19738/ *(results are not included in this scripts-only release; see manuscript archive)* contains the overall/sex-stratified DEG tables, GO-BP and targeted OXPHOS tables, summary, and report. It does not contain every intermediate RDS or every generated plot. |
| Path/run note | The script derives its dataset/cache paths from its location. The archive has not been configured as a portable clean-checkout run; this entry documents the current script contract, not a successful run from this repository. |

---

## 5. Key OXPHOS Results in v6.12

| Test | Value |
|------|-------|
| Preranked OXPHOS NES | −2.52 (downward; extracted from complete gseGO GO-BP run) |
| CAMERA fixed r = 0.01 (80 genes) | P = 1.77 × 10⁻⁹, within-method q = 5.30 × 10⁻⁹ |
| CAMERA estimated correlation (r = 0.178) | P = 0.0413, q = 0.124 across three microarrays |
| Directional FRY | P = 0.0439, q = 0.130 |
| Mean-z score base → +cell proxies | −0.218 → −0.171; adjusted nominal P = 0.0488, q = 0.195 |

**Interpretation**: This cohort showed the largest fixed-correlation OXPHOS enrichment, but after estimating inter-gene correlation from the data it remained only nominal and did not pass the three-microarray family correction. Direction of effect was consistently downward.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19738>
- Obtain deposited series matrix directly from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `GEOquery, clusterProfiler, data.table, limma, oligo, org.Hs.eg.db`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE19738__blood_microarray__differential_expression.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Spijker S, Van Zanten JS, De Jong S, et al. Stimulated gene expression profiles as a blood marker of major depressive disorder. *Biol Psychiatry*. 2010;68(2):179–186.
- Wu D, Smyth GK. Camera. *Nucleic Acids Res*. 2012;40:e133.

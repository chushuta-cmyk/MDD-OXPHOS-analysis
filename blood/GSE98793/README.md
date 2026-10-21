# GSE98793 — Whole-Blood Microarray Cohort (Discovery Screen + GAD Stratification)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE98793](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98793) |
| **Title** | Gene expression analysis in whole blood samples obtained from donors diagnosed with Major Depressive Disorder compared with healthy controls |
| **Also known as** | HiTDiP cohort (Leday et al. 2018) |
| **Species** | *Homo sapiens* |
| **Platform** | Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) |
| **Technology** | Microarray (deposited, GC-RMA normalised, log₂ scale) |
| **Reference** | Leday GGR, Vértes PE, Richardson S, et al. *Biol Psychiatry*. 2018;83(1):70–80. doi:10.1016/j.biopsych.2017.01.021 |

---

## 2. Sample Design

| Metric | Value |
|--------|-------|
| **Total analytic samples** | 192 |
| **MDD** | 128 (64 with comorbid generalised anxiety disorder / GAD; 64 without) |
| **Control** | 64 |
| **Tissue** | Whole blood |
| **Genes (post-QC, retained)** | 21,355 |
| **Covariates available** | Age, sex, batch, GAD status |
| **Note** | Two recorded batches, each containing 32 controls and 64 MDD participants |

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Primary exploratory blood GO-BP screen (cohort 1 of 2); also the GAD boundary-phenotype cohort |
| **Manuscript sections** | §2.1, §2.2, §2.3, §2.5; Table 1; Figure 1–2; Tables S2, S15–S18 |
| **Screening contribution** | Broad GO-BP returned 6,264 terms; 1,252 at BH q<0.05. OXPHOS NES = **−1.78, q = 0.00145**, which triggered the focused OXPHOS follow-up |
| **GAD stratification** | Re-modelled as Control / MDD-only / MDD+GAD (64 per group) with age+sex+batch; three contrasts tested for OXPHOS (Table S2) |

---

## 4. Processing in This Re-analysis

From Table S15:

- Loaded deposited processed expression matrix (median 5.713, consistent with log-scale input).
- Probe → symbol mapping via `hgu133plus2.db`; duplicate symbols collapsed to highest-mean probe.
- Gene-level design: `diagnosis + age + sex + batch`; robust/trend eBayes.
- Central OXPHOS pathway follow-up used the same matrix; `camera` called with `inter.gene.cor = 0.01` (fixed) and `= NA` (estimated); directional FRY as a self-contained comparator.
- Participant-level OXPHOS score: centered directional singscore, standardised within cohort; OLS with HC3 robust inference; covariate-nested models reported in Table S18.
- Cell-marker sensitivity: MCP-counter median-expression proxies added as covariates.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE98793__blood_microarray__differential_expression.R` |
| Required input | Local `GSE98793_series_matrix.txt`; if unavailable, the script falls back to GEOquery retrieval of GSE98793. `reference/MCPcounter/genes.txt` is also referenced for marker-based analyses. Neither file is included here. |
| Main outputs written by the script | `GSE98793_expression_cleaned.rds`, metadata CSV/RDS, all/significant DEG tables, three GAD contrast tables, MCP-counter/immune tables, GO/GSEA tables, volcano/QC figures, and report/session information. |
| Selected outputs included here | results/by_cohort/GSE98793/ *(results are not included in this scripts-only release; see manuscript archive)* contains DEG and GAD tables, metadata, immune/GSEA results, summary-related outputs, and volcano/report files. Intermediate expression RDS and all generated figures are not necessarily included. |
| Path/run note | The script uses absolute local input/output roots and an optional GEO network fallback. It is not portable from a clean checkout without path configuration. |

---

## 5. Key OXPHOS Results in v6.12

| Test | Value |
|------|-------|
| Preranked GO-BP OXPHOS NES (discovery) | −1.78, q = 0.00145 (downward) |
| CAMERA fixed r = 0.01 (136 genes) | P = 0.00249 (downward) |
| CAMERA estimated correlation (r = 0.137) | P = 0.295 (not significant) |
| Directional FRY | P = 0.358 |
| Rank-score β (age/sex/batch adjusted) | −0.126, HC3 P = 0.111, 95% CI −0.281 to 0.029 |
| Mean-z score base → +cell proxies | −0.040 → −0.020 (both null) |
| GAD contrasts (CAMERA estimated) | MDD+GAD vs Ctrl P = 0.115; MDD-only vs Ctrl P = 0.712; MDD+GAD vs MDD-only P = 0.348 |
| Three-cohort REML synthesis contribution | Included as screening cohort (k = 3) |

**Interpretation**: Under the default fixed-correlation assumption OXPHOS looked downward, but estimating inter-gene correlation from the same data removed the pathway-level signal. Participant-level scores were imprecise and batch dominated model R² (0.004 → 0.709 when batch added). No GAD contrast survived its correction family.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98793>
- Obtain the deposited series matrix / normalised expression directly from GEO under the original data-use terms.
- This repository does **not** redistribute participant-level expression.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `AnnotationDbi, GEOquery, clusterProfiler, dplyr, ggplot2, hgu133plus2.db, limma, matrixStats, org.Hs.eg.db, pheatmap, tidyr`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE98793__blood_microarray__differential_expression.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Leday GGR, Vértes PE, Richardson S, et al. Replicable and coupled changes in innate and adaptive immune gene expression in two case-control studies of blood microarrays in major depressive disorder. *Biol Psychiatry*. 2018;83(1):70–80.
- Ritchie ME et al. limma powers differential expression analyses. *Nucleic Acids Res*. 2015;43:e47.
- Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res*. 2012;40:e133.

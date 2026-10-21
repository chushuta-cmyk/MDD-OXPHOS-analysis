# GSE251778 — LeukoLOCK-Isolated Blood Leukocyte RNA-seq (Paired with GSE251780)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE251778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251778) |
| **Title** | Using multiomic integration to improve blood biomarkers of major depressive disorder: a case-control study (RNA-seq layer) |
| **Species** | *Homo sapiens* |
| **Platforms** | Illumina NovaSeq — **GPL24676** (female stratum) and **GPL20301** (male stratum); 3′-RNA-seq / bulk RNA-seq of LeukoLOCK-isolated peripheral blood leukocytes |
| **Technology** | RNA sequencing (raw counts, Ensembl-level) |
| **Reference** | Mokhtari A, Ibrahim EC, Gloaguen A, et al. *EBioMedicine*. 2025;113:105569. doi:10.1016/j.ebiom.2025.105569 |

---

## 2. Sample Design

| Stratum | Platform | MDD | Control | Total |
|---------|----------|-----|---------|-------|
| Female | GPL24676 | 50 | 55 | 105 |
| Male | GPL20301 | 30 | 34 | 64 |
| **Overall** | — | **80** | **89** | **169** |

- **Tissue**: LeukoLOCK-isolated peripheral blood leukocytes (not whole blood, not PBMC).
- **Covariates**: Sex (confounded with platform); diagnosis within joint sex/platform strata.
- **Paired methylation**: 147 of these 169 participants also have EPIC methylation in [GSE251780](../GSE251780/README.md).

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | RNA-seq follow-up cohort (separate family from the three microarrays); source of the paired RNA–promoter analysis |
| **Manuscript sections** | §2.1–§2.4, §2.6; Figure 1–4; Tables 1, 3, S6, S13–S15 |
| **Direction flag** | Only blood cohort with **upward** OXPHOS enrichment under broad GO-BP — opposite to the three microarrays. Sex and platform are completely confounded in this cohort. |

---

## 4. Processing in This Re-analysis

From Table S15:

- Ensembl version suffixes stripped; duplicate Ensembl IDs summed.
- Per-analysis filter: mean raw count > 10 within each stratum.
- Gene-level DESeq2 model: `sex + diagnosis` overall; `diagnosis` within the joint sex/platform strata.
- Ensembl-level rows retained for count models; symbol annotation kept separately for enrichment.
- Central pathway follow-up: separate `filterByExpr` + TMM + voom pipeline with diagnosis design; weighted CAMERA and FRY; unweighted sensitivity reported.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Per-cohort script | `GSE251778__blood_rnaseq__differential_expression.R` |
| Per-cohort inputs | `GSE251778_Mokhtari_rna_counts.csv.gz` and `GSE251778_metadata.csv`. The script expects these under `scripts/data/`; neither is included in this repository. |
| Per-cohort outputs | Under `scripts/results_fixed/`: `GSE251778_metadata_aligned.csv`, `GSE251778_DEG_{overall,female,male}.csv` and corresponding `*_adjP0.05.csv`, `GSE251778_GSEA_BP_*.csv`, `GSE251778_GO0006119_*.csv`, volcano PDFs, and `GSE251778_analysis_summary.csv` / report. |
| Weighted pathway sensitivity | `GSE251778__blood_rnaseq__weighted_camera_sensitivity.R` uses the count matrix and metadata; it writes `N2_GSE251778_weighted_camera_rebuild.tsv` and a session-info file to its configured cross-cohort results directory. |
| Selected outputs included here | results/by_cohort/GSE251778/ *(results are not included in this scripts-only release; see manuscript archive)* contains the DEG, GO-BP/OXPHOS, summary, metadata-alignment, and weighted CAMERA/FRY tables plus the overall volcano/report. The deposited count and metadata inputs are not included. |
| Path/run note | The primary script currently resolves its data/output roots relative to `scripts/`, while the weighted sensitivity script uses project-specific paths. These scripts are not yet portable from a clean GitHub checkout. |

---

## 5. Key OXPHOS Results in v6.12

| Test | Value |
|------|-------|
| Preranked GO-BP OXPHOS NES (overall) | **+1.99, q = 2.52 × 10⁻⁴** (upward) |
| Preranked OXPHOS NES, male/GPL20301 | +2.56, q = 8.14 × 10⁻¹² (upward) |
| Preranked OXPHOS NES, female/GPL24676 | −1.29, q = 0.508 (downward, null) |
| Weighted CAMERA fixed r = 0.01 (136 genes) | P = 0.00651 (upward) |
| Weighted CAMERA estimated correlation (r = 0.266) | P = 0.492 (null) |
| Weighted directional FRY | P = 0.182 |
| Rank-score β (four-cohort sensitivity) | Included; point estimate positive |
| Mean-z score base → +cell proxies | +0.082 → +0.109 (both null) |
| Paired 147-subset RNA score, before → after DNAm cell adjustment | +0.080 (P = 0.580) → +0.127 (P = 0.427) |

**Interpretation**: The upward OXPHOS signal was almost entirely attributable to the male/GPL20301 stratum under fixed correlation; estimating inter-gene correlation (r = 0.266) collapsed it. Sex and platform cannot be disentangled here, so this cohort is not an independent replication of the microarray downward trend.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251778>
- RNA raw counts available from GEO. This repository does not redistribute participant-level data.
- Methylation layer for the same participants: see [GSE251780](../GSE251780/README.md).

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `AnnotationDbi, DESeq2, data.table, edgeR, ggplot2, limma, org.Hs.eg.db`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE251778__blood_rnaseq__differential_expression.R
Rscript ./GSE251778__blood_rnaseq__weighted_camera_sensitivity.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Mokhtari A, Ibrahim EC, Gloaguen A, et al. Using multiomic integration to improve blood biomarkers of major depressive disorder: a case-control study. *EBioMedicine*. 2025;113:105569.
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol*. 2014;15:550.
- Law CW et al. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol*. 2014;15:R29.

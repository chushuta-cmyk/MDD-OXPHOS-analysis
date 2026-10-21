# GSE80655 — Postmortem Brain RNA-seq (Three Regions, Psychiatric Cohort)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE80655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655) |
| **Title** | Post-mortem molecular profiling of three psychiatric disorders |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina HiSeq 2000; NCBI GRCh38.p13 raw counts |
| **Technology** | Bulk RNA sequencing (postmortem brain) |
| **Reference** | Ramaker RC, Bowling KM, Lasseigne BN, et al. *Genome Med*. 2017;9:72. doi:10.1186/s13073-017-0458-5 |

---

## 2. Sample Design

| Region | Analytic N (MDD/control) |
|--------|-------------------------:|
| Anterior cingulate gyrus (AnCg) | 48 |
| Nucleus accumbens (nAcc) | 44 |
| Dorsolateral prefrontal cortex (DLPFC) | 47 |

- Original cohort also includes schizophrenia and bipolar disorder; only MDD vs control retained.
- **Genes retained**: ~35,000 per region (Entrez → symbol, largest summed count; count ≥1 in ≥10 samples).

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Three of seven regional brain OXPHOS comparisons (RNA-seq); AnCg carries the strongest estimated-correlation brain signal |
| **Manuscript sections** | §2.2–§2.3; Figure 3; Table 2; Tables S2, S15–S16 |

---

## 4. Processing in This Re-analysis

From Table S15:

- Regional raw counts; Entrez → symbol, largest summed count; count ≥1 in ≥10 samples.
- Gene-level DESeq2 model: `diagnosis`; shrunken effects kept separate from output statistics.
- **Central pathway follow-up**: separate `filterByExpr` + TMM + voom pipeline with diagnosis design; weighted precision weights retained.
- OXPHOS pathway: `camera` fixed r = 0.01 and estimated r; directional FRY.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE80655__brain_rnaseq__differential_expression.R` |
| Required input | `GSE80655_raw_counts_GRCh38.p13_NCBI.tsv` (NCBI count matrix). Optional local phenotype inputs are `GSE80655_series_matrix.txt` and/or `GSE80655_metadata.csv` or `.rds`; GEOquery is a fallback for phenotype retrieval. These source files are not included here. |
| Main outputs written by the script | Per-region `GSE80655_all_DE_results_<region>.csv`, merged `GSE80655_MDD_vs_Control_BrainRegion_DEGs.csv`, regional GO-BP/KEGG/Reactome tables, QC/PCA/correlation and volcano plots, MOFA expression/metadata/feature-annotation RDS files, analysis summary, and session information. |
| Selected outputs included here | results/by_cohort/GSE80655/ *(results are not included in this scripts-only release; see manuscript archive)* contains per-region and merged DEG tables, analysis summary, volcano PDFs, and later SLC6A2 diagnostic/sensitivity tables. Raw counts and MOFA intermediates are not included. |
| Path/run note | The script writes to the original workstation's `transcriptomics/brain/GSE80655/results_fixed` path and has a legacy fallback metadata path. The old standalone README's four-step workflow is not the script structure present in this archive; the linked single integrated script is the current archived code entry. |

---

## 5. Key OXPHOS Results in v6.12

| Region | Fixed CAMERA P | Estimated CAMERA P (r) | FRY P | FRY q (7 regions) |
|--------|---------------:|-----------------------:|------:|------------------:|
| AnCg | 2.16 × 10⁻¹⁴ | **0.0642 (r = 0.275)** | 0.0267 | 0.187 |
| nAcc | — | 0.219 | 0.183 | 0.641 |
| DLPFC | — | 0.894 | 0.950 | 0.950 |

- Direction in AnCg was **upward** under both settings.
- Across all 11 brain+blood comparisons, estimated inter-gene correlations ranged 0.063–0.343; only one comparison remained nominal under estimated correlation, and none passed its available family correction.

**Interpretation**: The AnCg upward OXPHOS signal was highly sensitive to the inter-gene-correlation assumption (P from 2 × 10⁻¹⁴ to 0.064); direction did not flip.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655>
- Obtain deposited raw count matrix from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `AnnotationDbi, DESeq2, GEOquery, RColorBrewer, ReactomePA, apeglm, clusterProfiler, dplyr, enrichplot, ggplot2, org.Hs.eg.db, pheatmap, tibble`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE80655__brain_rnaseq__differential_expression.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Ramaker RC, Bowling KM, Lasseigne BN, et al. Post-mortem molecular profiling of three psychiatric disorders. *Genome Med*. 2017;9:72.
- Love MI, Huber W, Anders S. DESeq2. *Genome Biol*. 2014;15:550.
- Law CW et al. voom. *Genome Biol*. 2014;15:R29.

# GSE53987 — Postmortem Brain Microarray (Three Regions, Tetrad Design)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE53987](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987) |
| **Title** | Postmortem transcriptional profiling reveals widespread increase in inflammation in schizophrenia: a comparison of prefrontal cortex, striatum, and hippocampus among matched tetrads of controls with subjects diagnosed with schizophrenia, bipolar or major depressive disorder |
| **Species** | *Homo sapiens* |
| **Platform** | Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) |
| **Technology** | Microarray (postmortem brain) |
| **Reference** | Lanz TA, Reinhart V, Sheehan MJ, et al. *Transl Psychiatry*. 2019;9:151. doi:10.1038/s41398-019-0492-8 |

---

## 2. Sample Design

| Region | Analytic N | Design |
|--------|-----------:|--------|
| Hippocampus (HIP) | 32 | MDD vs control, within matched tetrads |
| Prefrontal cortex BA46 (DLPFC) | 34 | MDD vs control, within matched tetrads |
| Associative striatum (STR) | 32 | MDD vs control, within matched tetrads |

- Original tetrad design also included schizophrenia and bipolar disorder; only MDD and control samples retained.
- **Genes (post-QC)**: 21,355 per region.
- **Covariates**: Diagnosis, age, postmortem interval (PMI), pH, matched tetrad block.

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Three of seven regional brain OXPHOS comparisons (microarray) |
| **Manuscript sections** | §2.2–§2.3; Figure 3; Table 2; Tables S2, S15–S16 |
| **Note** | Original matrices were already log-scale; no added quantile normalisation in v6.12. |

---

## 4. Processing in This Re-analysis

From Table S15:

- Loaded deposited matrix; `log2(x+1)` only if needed; expression > 3 in ≥10% samples.
- Probe → symbol via `hgu133plus2.db`; duplicate symbols collapsed to highest-mean probe.
- Region-stratified design: `diagnosis + age + PMI + pH + tetrad`; ordinary eBayes.
- OXPHOS pathway follow-up: `camera` fixed r = 0.01 and estimated r; directional FRY.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE53987__brain_microarray__differential_expression.R` |
| Required input | `GSE53987_series_matrix.txt` (deposited expression and sample annotations); the script checks local relative paths and then workstation-specific fallback paths. The matrix is not included here. |
| Main outputs written by the script | Per-region `GSE53987_all_DE_results_<region>.csv`, merged `GSE53987_MDD_vs_Control_BrainRegion_DEGs.csv`, region expression/clinical CSVs, MOFA expression/metadata RDS files, regional volcano PDFs, `GSE53987_analysis_summary.csv`, and session information. |
| Selected outputs included here | results/by_cohort/GSE53987/ *(results are not included in this scripts-only release; see manuscript archive)* contains the per-region and merged DEG tables, clinical tables, summary, and volcano PDFs. The script's MOFA RDS outputs are not included. |
| Path/run note | Output is currently hard-coded to the original workstation's `transcriptomics/brain/GSE53987/results_fixed` directory. The script also uses an expression-range rule to decide whether to apply `log2(x+1)`; the README does not infer the deposited matrix scale independently. |

---

## 5. Key OXPHOS Results in v6.12

| Region | CAMERA estimated P | FRY P | Direction |
|--------|-------------------:|------:|-----------|
| Hippocampus | 0.406 | 0.565 | — |
| BA46 / DLPFC | 0.829 | 0.902 | — |
| Associative striatum | 0.676 | 0.494 | — |

- Preranked hippocampal NES was the most downward among brain regions (NES ≈ −2.11), but did not survive estimated-correlation correction.
- None of the three regions contributed a nominal estimated-correlation result.

**Interpretation**: No regional OXPHOS signal survived correlation-aware testing in this cohort.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987>
- Obtain deposited series matrix directly from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `AnnotationDbi, edgeR, hgu133plus2.db, limma, stringr`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE53987__brain_microarray__differential_expression.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Lanz TA, Reinhart V, Sheehan MJ, et al. Postmortem transcriptional profiling reveals widespread increase in inflammation in schizophrenia. *Transl Psychiatry*. 2019;9:151.
- Wu D, Smyth GK. Camera. *Nucleic Acids Res*. 2012;40:e133.

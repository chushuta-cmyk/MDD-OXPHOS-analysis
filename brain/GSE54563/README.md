# GSE54563 — Postmortem Anterior Cingulate Cortex Microarray (Matched Pairs)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE54563](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54563) |
| **Title** | Gene expression in anterior cingulate cortex in major depressive disorder |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina HumanHT-12 V3.0 expression beadchip (GPL6947) |
| **Technology** | Microarray (postmortem brain) |
| **Reference** | Chang LC, Jamain S, Lin CW, et al. *PLoS One*. 2014;9(3):e90980. doi:10.1371/journal.pone.0090980 |

---

## 2. Sample Design

| Metric | Value |
|--------|-------|
| **Total samples** | 50 |
| **MDD** | 25 |
| **Control** | 25 |
| **Design** | Matched pairs |
| **Brain region** | Anterior cingulate cortex (ACC) |
| **Probes (raw)** | 48,803 |
| **Genes (post-QC, retained)** | 19,478 |
| **Covariates** | Diagnosis, matched pair |

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | One of seven regional brain OXPHOS comparisons (microarray, matched-pair design) |
| **Manuscript sections** | §2.2–§2.3; Figure 3; Table 2; Tables S2, S15–S16 |

---

## 4. Processing in This Re-analysis

From Table S15:

- Loaded deposited matrix; `log2(x+1)`; expression > 3 in ≥10% samples.
- Probe-level BH correction across all 48,803 probes **before** annotation; representative probe per gene selected by highest average expression.
- Design: `diagnosis + pair`; ordinary eBayes.
- OXPHOS pathway: `camera` fixed r = 0.01 and estimated r; directional FRY.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE54563__brain_microarray__differential_expression.R` |
| Required input | `GSE54563_series_matrix.txt` (contains the expression matrix and matched-pair sample annotations); the script searches local paths including `data/`. The input is not included here. |
| Main outputs written by the script | `GSE54563_All_DEG_Results.csv`, `GSE54563_Significant_DEGs_adjP0.05.csv`, GO/KEGG/Reactome enrichment tables and plots, volcano TIFF/PDF, analysis summary, and MOFA expression/metadata files. |
| Selected outputs included here | results/by_cohort/GSE54563/ *(results are not included in this scripts-only release; see manuscript archive)* contains the all/significant DEG tables, analysis summary, and volcano PDF. The MOFA intermediate files and enrichment outputs are not included in this selected archive. |
| Path/run note | Output is hard-coded to the original workstation's `transcriptomics/brain/GSE54563/results_fixed` directory; the repository copy has not been path-normalized or test-run. |

---

## 5. Key OXPHOS Results in v6.12

| Test | Value |
|------|-------|
| CAMERA estimated correlation | P = 0.64 |
| Directional FRY | P = 0.72 |

No nominal OXPHOS association in ACC for this cohort.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54563>
- Obtain deposited series matrix directly from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `AnnotationDbi, ReactomePA, clusterProfiler, data.table, dplyr, enrichplot, ggplot2, illuminaHumanv3.db, limma, matrixStats, org.Hs.eg.db`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE54563__brain_microarray__differential_expression.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Chang LC, Jamain S, Lin CW, et al. A conserved BDNF, glutamate- and GABA-enriched gene module related to human depression identified by coexpression meta-analysis and DNA variant genome-wide association studies. *PLoS One*. 2014;9:e90980.
- Wu D, Smyth GK. Camera. *Nucleic Acids Res*. 2012;40:e133.

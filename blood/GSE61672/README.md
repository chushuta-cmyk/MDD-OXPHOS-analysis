# GSE61672 — Community Blood Microarray for Generalised Anxiety Symptoms

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE61672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61672) |
| **Title** | Blood gene expression profiles suggest altered immune function associated with symptoms of generalized anxiety disorder |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina HumanHT-12 V4.0 expression beadchip (GPL10558) |
| **Technology** | Microarray (whole blood, community sample) |
| **Reference** | Wingo AP, Gibson G. *Brain Behav Immun*. 2015;43:184–191. doi:10.1016/j.bbi.2014.09.016 |

---

## 2. Sample Design

| Metric | Value |
|--------|-------|
| **GEO total samples** | 546 community participants |
| **Analytic in v6.12** | 531 (continuous GAD-7 model) and 336 (anxiety case–control model), complete cases per model |
| **Tissue** | Whole blood |
| **Covariates** | Age, BMI, ethnicity, batch, RNA integrity number (RIN) |

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Boundary phenotype: continuous GAD-7 and anxiety case–control OXPHOS tests |
| **Manuscript sections** | §2.5; Table S2; Tables S15–S16 |
| **Note** | Two distinct models (continuous and binary) are corrected separately; no pooled anxiety replication. |

---

## 4. Processing in This Re-analysis

From Table S15:

- Loaded deposited expression; `illuminaHumanv4.db` annotation; duplicate symbols collapsed to highest-mean probe.
- Complete cases per model.
- Design: `GAD-7 score or case status + age + BMI + ethnicity + batch + RIN`; robust/trend eBayes.
- OXPHOS pathway: `camera` fixed r = 0.01 and estimated r; directional FRY; correction within the two anxiety comparisons.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE61672__blood_microarray__anxiety_boundary.R` |
| Required input | `GSE61672_series_matrix.txt.gz`, expected under `analysis_framework/GSE61672/data/`; the script also sources `analysis_framework/_libpaths.R`. Neither is included here. |
| Main outputs written by the script | `GSE61672_metadata_parsed.csv`, `GSE61672_gene_expression.rds`, `GSE61672_metadata.rds`; `GSE61672_DE_{continuous_GAD7,case_vs_control}.tsv`; `GSE61672_OXPHOS_{CAMERA,fgsea,score_models}.tsv`; `GSE61672_sample_scores.tsv`; `GSE61672_analysis_summary.txt`; session info; scatter/case-control/forest/volcano/PCA PNG and PDF figures. |
| Selected outputs included here | results/by_cohort/GSE61672/ *(results are not included in this scripts-only release; see manuscript archive)* contains the two DE tables, pathway/score tables, sample scores, and summary. Input matrix and figures are not included there. |
| Path/run note | The script expects an `analysis_framework/` tree and a configured local library-path helper; a clean repository checkout does not currently provide those dependencies. |

---

## 5. Key OXPHOS Results in v6.12

| Contrast | CAMERA estimated P | FRY P | FRY q (2 comparisons) |
|----------|-------------------:|------:|----------------------:|
| GAD-7 continuous (N = 531) | 0.768 | 0.554 | 0.554 |
| Anxiety case vs control (N = 336) | 0.546 | 0.307 | 0.554 |

No anxiety contrast approached the OXPHOS correction family.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61672>
- Obtain deposited series matrix directly from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `AnnotationDbi, Biobase, GEOquery, data.table, fgsea, ggplot2, ggrepel, illuminaHumanv4.db, limma, lmtest, org.Hs.eg.db, sandwich, singscore`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE61672__blood_microarray__anxiety_boundary.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Wingo AP, Gibson G. Blood gene expression profiles suggest altered immune function associated with symptoms of generalized anxiety disorder. *Brain Behav Immun*. 2015;43:184–191.
- Wu D, Smyth GK. Camera. *Nucleic Acids Res*. 2012;40:e133.

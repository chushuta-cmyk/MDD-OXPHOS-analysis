# GSE228702 — Blood 3′ RNA-seq for Social Anxiety and Early-Life Adversity

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE228702](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228702) |
| **Title** | Blood transcriptome analysis suggests an indirect molecular association of early life adversities (ELA) and adult social anxiety disorder (SAD) by immune-related signal transduction |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina NovaSeq 6000 (GPL24247); 3′-RNA-seq (Lexogen) |
| **Technology** | Bulk RNA sequencing (whole blood) |
| **Reference** | Edelmann S, Wiegand A, Hentrich T, et al. *Front Psychiatry*. 2023;14:1125553. doi:10.3389/fpsyt.2023.1125553 |

---

## 2. Sample Design

| Metric | Value |
|--------|-------|
| **Analytic samples in v6.12** | 160 |
| **Design** | Factorial: social anxiety disorder (SAD) × early-life adversity (ELA), four groups |
| **Tissue** | Whole blood |
| **Covariates** | Age, gender |
| **Input matrix** | Author-deposited cell-adjusted expression matrix (not raw counts) |

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Boundary phenotype: SAD marginal, ELA marginal, and SAD×ELA interaction OXPHOS tests |
| **Manuscript sections** | §2.5; Figure 1; Table S2; Tables S15–S16 |

---

## 4. Processing in This Re-analysis

From Table S15:

- Used the deposited cell-adjusted expression matrix directly.
- Largest-variance symbol retained for duplicate probes.
- Within-gene rank inverse-normal transformation applied for DE and pathway tests.
- Design: `0 + four groups + age + gender`; SAD and ELA marginal contrasts plus their interaction.
- Gene-level limma, followed by CAMERA/FRY with the corresponding design.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE228702__blood_rnaseq__sad_ela_boundary.R` |
| Required inputs | `GSE228702_adjusted_expression_Cellcounts_granulatorAbis0_nnls.csv.gz` and `GSE228702_series_matrix.txt.gz`, expected under `analysis_framework/GSE228702/data/`. Neither is included here. The script also sources `analysis_framework/_libpaths.R`. |
| Main outputs written by the script | `GSE228702_metadata_parsed.csv`; `GSE228702_DE_{SAD_marginal,ELA_marginal,SADxELA_interaction}.tsv`; `GSE228702_OXPHOS_{CAMERA,fgsea,score_models}.tsv`; `GSE228702_sample_scores.tsv`; `GSE228702_analysis_summary.txt`; session info; group/forest/volcano/PCA PNG and PDF figures. |
| Selected outputs included here | results/by_cohort/GSE228702/ *(results are not included in this scripts-only release; see manuscript archive)* contains the three DE contrasts, pathway/score tables, sample scores, and summary. The adjusted matrix and figures are not included in that directory. |
| Path/run note | The script expects to be run with a project root containing `analysis_framework/`; that structure is not present in this repository package. |

---

## 5. Key OXPHOS Results in v6.12

| Contrast | CAMERA estimated P | FRY P | FRY q (3 contrasts) |
|----------|-------------------:|------:|--------------------:|
| SAD main contrast | 0.583 | 0.775 | 0.775 |
| Adversity main contrast | 0.115 | 0.626 | 0.775 |
| SAD × adversity interaction | 0.968 | 0.364 | 0.775 |

No social-anxiety or adversity contrast approached OXPHOS correction.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228702>
- Obtain deposited expression matrix from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `AnnotationDbi, Biobase, GEOquery, data.table, fgsea, ggplot2, ggrepel, limma, org.Hs.eg.db, sandwich, singscore`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE228702__blood_rnaseq__sad_ela_boundary.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Edelmann S, Wiegand A, Hentrich T, et al. Blood transcriptome analysis suggests an indirect molecular association of early life adversities and adult social anxiety disorder by immune-related signal transduction. *Front Psychiatry*. 2023;14:1125553.
- Wu D, Smyth GK. Camera. *Nucleic Acids Res*. 2012;40:e133.

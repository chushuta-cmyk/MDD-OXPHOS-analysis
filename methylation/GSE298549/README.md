# GSE298549 — Whole-Blood 450K Methylation by Weighted Life Events (Exploratory)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE298549](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298549) |
| **Title** | Effects of stressful life-events on DNA methylation in panic disorder and major depressive disorder (MDD cohort, wLE analysis) |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina HumanMethylation450 BeadChip |
| **Technology** | DNA methylation array (whole blood) |
| **Reference** | Czamara D, Neufang A, Dieterle R, et al. *Clin Epigenetics*. 2022;14:55. doi:10.1186/s13148-022-01274-y |

---

## 2. Sample Design

| Metric | Value |
|--------|-------|
| **Analytic samples** | 102 MDD patients; 74 with complete weighted-life-event (wLE) records |
| **Retained probes** | 320,481 (autosomal, finite in ≥95% of samples) |
| **Tissue** | Whole blood |
| **Phenotype** | Standardised weighted life-event score (continuous), not case–control MDD status |

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Exploratory external methylation model (Supplementary Text S8); not pooled with the paired GSE251780 clinical-MDD promoter analysis |
| **Manuscript sections** | Supplementary Text S8; Table S28; Supplementary Data S16 |

---

## 4. Processing in This Re-analysis

From Supplementary Text S8:

- Processed blood M-value matrix.
- Autosomal probes finite in ≥95% of samples → 320,481 tests.
- Design: `standardised weighted life events + standardised age + sex + array position`; limma eBayes with robust and trend settings.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE298549__blood_methylation__weighted_life_events.R` |
| Required inputs | `GSE298549_mval_MDD_cases.txt.gz` and `GSE298549_metadata.csv` under the script's configured `analysis_framework/GSE201016/data/` directory. The directory name is a legacy path label; the input accession is GSE298549. These files are not included here. The script also sources `analysis_framework/_libpaths.R`. |
| Main outputs written by the script | `GSE298549_wLE_EWAS_{all,FDR05}.tsv`, `GSE298549_wLE_EWAS_sparse_slide_sensitivity.tsv`, `GSE298549_reported_CpG_replication.tsv`, `GSE298549_OXPHOS_promoter_probes.txt`, `GSE298549_OXPHOS_gsameth.tsv`, `GSE298549_OXPHOS_score_model.tsv`, `GSE298549_OXPHOS_sample_scores.tsv`, optional `GSE298549_wLE_DMRcate.tsv`, figures, summary, and session information. |
| Selected outputs included here | results/by_cohort/GSE298549/ *(results are not included in this scripts-only release; see manuscript archive)* contains the FDR EWAS, reported-CpG lookup, gsameth and score tables, and summary. The full EWAS and figures are not included in this selected archive. |
| Path/run note | The script writes under `analysis_framework/GSE201016/` despite analyzing GSE298549; this legacy path has not been renamed in the code. |

---

## 5. Key OXPHOS Results in v6.12

- Reported only as an exploratory external methylation model; no effects pooled with the paired clinical-MDD promoter tests.
- The saved OXPHOS score CI is a normal interval; selected single-CpG intervals are reconstructed moderated-t intervals (not OLS or Δβ intervals).

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298549>
- Obtain deposited methylation data from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `IlluminaHumanMethylation450kanno.ilmn12.hg19, data.table, ggplot2, ggrepel, limma, lmtest, minfi, missMethyl, org.Hs.eg.db, sandwich`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE298549__blood_methylation__weighted_life_events.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Czamara D, Neufang A, Dieterle R, et al. Effects of stressful life-events on DNA methylation in panic disorder and major depressive disorder. *Clin Epigenetics*. 2022;14:55.

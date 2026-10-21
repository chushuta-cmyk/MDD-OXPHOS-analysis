# GSE251780 — Paired Blood DNA Methylation (Illumina EPIC, Subset of GSE251778)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE251780](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251780) |
| **Title** | Using multiomic integration to improve blood biomarkers of major depressive disorder: a case-control study (DNA methylation layer) |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina HumanMethylationEPIC BeadChip |
| **Technology** | DNA methylation array (processed beta values) |
| **Reference** | Mokhtari A, Ibrahim EC, Gloaguen A, et al. *EBioMedicine*. 2025;113:105569. doi:10.1016/j.ebiom.2025.105569 |

---

## 2. Sample Design

| Metric | Value |
|--------|-------|
| **Total analytic participants** | 147 (all nested within the 169 of [GSE251778](../GSE251778/README.md)) |
| **Tissue** | Peripheral blood leukocytes (same LeukoLOCK isolates as RNA) |
| **Retained CpGs** | 723,975 (after detection-P, autosomal, common CpG/SBE-SNP filtering) |
| **Covariates** | Sex; six EpiDISH RPC cell estimates (B, NK, CD4T, CD8T, Mono, Eosino; neutrophils omitted in the cell-adjusted model) |

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Paired RNA–DNAm cohort; sole source of the EWAS and promoter-disease association tests |
| **Manuscript sections** | §2.6; Figure 5; Tables S6, S14, S19–S21, S27; Supplementary Data S14–S15 |
| **Not a separate RNA cohort** | The 147 paired participants are counted inside GSE251778, not as an additional independent cohort. |

---

## 4. Processing in This Re-analysis

From §2.6 and Table S27:

- Processed EPIC beta values filtered on detection P, autosomal annotation, and common CpG/SBE SNP criteria → 723,975 CpGs.
- M-value regression with `limma`:
  - **Sex-adjusted model**: diagnosis + sex.
  - **Cell-adjusted model**: diagnosis + sex + six EpiDISH RPC cell estimates (neutrophils omitted to avoid compositional collinearity).
- Promoter scores: TSS1500 / TSS200 / first-exon annotations; average M values within genes, standardise genes, average covered genes, standardise set score.
- Seven gene sets (GO OXPHOS, nuclear MitoCarta respiratory subunits, complexes I–V) × two models = **14 promoter-score tests**.
- Probe-bias-aware set test: `gsameth` / missMethyl on the retained EPIC universe, GO OXPHOS Entrez collection, raw-CpG cutoffs 0.001 / 0.005 / 0.01.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Paired EWAS/RNA-DNAm script | `GSE251780__blood_methylation__paired_ewas_rna_dnam.R` |
| Cell-adjusted sensitivity script | `GSE251780__blood_methylation__cell_adjusted_sensitivity.R` |
| Primary inputs | `GSE251780_Mokhtari_MatrixProcessed_AverageAndPvalue.txt.gz`, `GSE251780_metadata.csv`, the existing `GSE251780_MDD_DMP_all.tsv` and `GSE251780_paired_sample_scores.tsv`, plus GSE251778 counts `GSE251778_Mokhtari_rna_counts.csv.gz`. The script also consumes the existing cross-cohort score file `blood_OXPHOS_sample_scores.tsv`. These source inputs are not distributed here. |
| Main outputs written by paired script | `GSE251780_MDD_DMP_all.tsv`, `GSE251780_MDD_DMP_FDR05.tsv`, `GSE251780_detectionP_sample_QC.tsv`, `GSE251780_MDD_OXPHOS_gsameth.tsv`, `GSE251780_OXPHOS_DNAm_score_model.tsv`, `GSE251780_paired_pathway_association.tsv`, `GSE251780_paired_sample_scores.tsv`, `GSE251780_paired_OXPHOS_gene_associations.tsv`, `GSE251780_paired_gene_set_permutation.tsv`, figures, summary, and session information. |
| Main outputs written by cell-adjustment script | `GSE251780_EpiDISH_cell_fractions.tsv`, `GSE251780_MDD_DMP_cell_adjusted_all.tsv`, `GSE251780_MDD_DMP_cell_adjusted_FDR05.tsv`, `GSE251780_OXPHOS_DNAm_score_cell_adjusted.tsv`, `GSE251780_paired_pathway_cell_adjusted.tsv`, `GSE251780_paired_pathway_cell_adjusted_by_sex.tsv`, `GSE251780_DNAJC15_cell_adjusted.tsv`, `GSE251780_EpiDISH_group_differences.tsv`, `GSE251780_cell_adjusted_summary.txt`, and a paired RNA–DNAm figure. |
| Promoter follow-up scripts | Historical Plan A: `GSE251780__blood_methylation__promoter_planA_historical.R`; v3.22 targeted recalculation: `GSE251780__blood_methylation__promoter_v322_recalculation.R`; saved-table assembly only: `GSE251780__blood_methylation__promoter_v322_table_assembly.R`; Figure 5 export: `GSE251780__blood_methylation__figure5_promoter_v322.R`. The v3.22 table-assembly script combines saved inputs and does not fit models. |
| Promoter-specific inputs | The historical Plan A workflow uses `Human.MitoCarta3.0.xls`, existing EWAS/probe data, metadata, and cell fractions. The v3.22 recalculation uses `prepared_mapping.rds`, `target_promoter_beta.tsv`, `GSE251780_metadata.csv`, and `GSE251780_EpiDISH_cell_fractions.tsv`; table assembly additionally reads saved Plan A and v3.22 outputs. Several are under external `outputs/secondary_planA_2026-10-05/` or `analysis_framework/GSE251780/` paths and are not included here. |
| Selected outputs included here | results/by_cohort/GSE251780/ *(results are not included in this scripts-only release; see manuscript archive)* contains the listed paired EWAS, score, gsameth, cell-fraction, paired-association, promoter, and audit tables. The full processed EPIC/RNA inputs and intermediate mapping objects are not included. Figure 5 source data and rendered outputs are also present in the repository-level `figure_source_data/` and `figures/` directories. |
| Path/run note | These scripts use a mix of `analysis_framework/`, `transcriptomics/`, `outputs/secondary_planA_2026-10-05/`, and absolute workstation paths. They document distinct historical/primary/sensitivity steps and are not a one-command portable pipeline in the current archive. |

---

## 5. Key OXPHOS Results in v6.12

| Analysis | Result |
|----------|--------|
| Full EWAS, sex-adjusted | 0 CpGs at BH q < 0.05 |
| Full EWAS, cell-adjusted | 0 CpGs at BH q < 0.05 |
| GO OXPHOS promoter score, sex-adjusted β | 0.036 (95% CI −0.299 to 0.371; P = 0.833) |
| GO OXPHOS promoter score, cell-adjusted β | 0.061 (95% CI −0.285 to 0.406; P = 0.729) |
| Nuclear MitoCarta subunits, sex → cell adjusted | 0.025 → 0.049 (both null) |
| Complexes I–V, cell-adjusted range | 0.0005 (CIV) to 0.076 (CI); all 95% CI cross zero |
| All 14 promoter tests | BH q = 0.998 for every test |
| gsameth P at CpG cutoff 0.001 / 0.005 / 0.01 | 1.000 / 0.921 / 0.873 |
| MDD cell-fraction differences | Neutrophils +0.0239 (P = 0.190); NK −0.0090 (P = 0.175); none survived q = 0.622 |

**Interpretation**: No corrected OXPHOS promoter or CpG-level disease association in this paired blood cohort. Effects are on standardised gene-averaged promoter M-value scores, not on Δβ of individual CpGs.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251780>
- EPIC beta values available from GEO under the original data-use terms. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `IlluminaHumanMethylationEPICanno.ilm10b4.hg19, data.table, edgeR, ggplot2, ggrepel, limma, lmtest, minfi, missMethyl, org.Hs.eg.db, readxl, sandwich, utils`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE251780__blood_methylation__cell_adjusted_sensitivity.R
Rscript ./GSE251780__blood_methylation__paired_ewas_rna_dnam.R
Rscript ./GSE251780__blood_methylation__promoter_planA_historical.R
Rscript ./GSE251780__blood_methylation__promoter_v322_recalculation.R
Rscript ./GSE251780__blood_methylation__promoter_v322_table_assembly.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Mokhtari A, Ibrahim EC, Gloaguen A, et al. Using multiomic integration to improve blood biomarkers of major depressive disorder. *EBioMedicine*. 2025;113:105569.
- Teschendorff AE et al. A comparison of reference-based algorithms for correcting cell-type heterogeneity in epigenome-wide association studies. *BMC Bioinformatics*. 2017;18:105.
- Maksimovic J, Oshlack A, Phipson B. Gene set enrichment analysis for genome-wide DNA methylation data. *Genome Biol*. 2021;22:173.

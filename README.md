# MDD–OXPHOS Multi-Omics Re-analysis — Analysis Scripts

Analysis pipeline for *"Reassessing oxidative phosphorylation associations in depression using public transcriptomic, methylomic and genetic data"* (retrospective exploratory re-analysis). This repository contains **only the analysis scripts** (31 R / Python pipelines) and this README; raw data and result tables are **not** distributed.

All 31 pipelines were executed end-to-end on the author's workstation on 2026-10-21 (round-1 run-through: only run-blocking path/compatibility issues were fixed — no statistical model, threshold, gene set or result value was altered). For this public snapshot, every machine-specific absolute path in the scripts has been replaced by **placeholder roots** (see [Path configuration](#path-configuration)); after that normalisation a full syntax check of all 31 scripts was re-run. The scripts are archived analysis code: to re-run them you must edit the placeholder paths and supply the inputs described below.

## Quick start

1. **Install R (≥ 4.5)** and the packages listed in [Environment](#environment); install the Python packages for the `genetics/` and `lookups/` scripts.
2. **Obtain the input files.** Each pipeline reads from one of the roots described below: `DATA_ROOT` (raw GEO / supplementary downloads — see [Inputs / Outputs per pipeline](#inputs--outputs-per-pipeline)), `WORK_ROOT` and `ARCHIVE_ROOT` (author-produced intermediate result tables from earlier runs — these are **not** on GEO and not shipped here).
3. **Edit the placeholder paths.** Every script carries a short "Path configuration" block near the top defining `DATA_ROOT`, `WORK_ROOT`, `ARCHIVE_ROOT` and `OUTPUT_ROOT` as `<PLACEHOLDER>` strings. Replace them with the actual directories on your machine (see [Path configuration](#path-configuration)).
4. **Run from the repository root:**

   ```bash
   Rscript blood/GSE98793/GSE98793__blood_microarray__differential_expression.R
   python3 genetics/MULTI__genetics__fuma_mdd_anxiety.py
   ```

   Results are written under `OUTPUT_ROOT` (default `./output/`), with per-pipeline subpaths listed in the table below.

## Path configuration

Each script defines four placeholder roots (edit them before running):

| Variable | Meaning | Placeholder example |
|---|---|---|
| `DATA_ROOT` | Directory holding the raw GEO / supplementary input files (series matrices, count tables, manifests, normalised beta matrices, metadata CSVs). Expected file names per pipeline are listed in the table below. | `<DATA_ROOT>` |
| `WORK_ROOT` | Author's analysis workspace holding intermediate result tables consumed by the `MULTI__` cross-cohort, `lookups/` and `methylation/` pipelines. Not redistributed — supply from your own previous runs. | `<WORK_ROOT>` |
| `ARCHIVE_ROOT` | Author's assembled manuscript output archive consumed by `figures/`, `genetics/` and some `lookups/`/`methylation/` pipelines. Not redistributed. | `<ARCHIVE_ROOT>` |
| `OUTPUT_ROOT` | Where this run writes results. Default `./output` (relative to the repository root). | `./output` |

All paths below are **logical** — e.g. `DATA_ROOT/GSE53987/GSE53987_series_matrix.txt` means the file `GSE53987_series_matrix.txt` under the directory you set as `DATA_ROOT` (scripts search several candidate sub-locations and print the path actually found).

## Naming convention

Scripts follow `GSE<accession>__<tissue>_<modality>__<purpose>.<ext>` for single-accession analyses and `MULTI__<scope>__<purpose>.<ext>` for cross-cohort analyses. The old-name → new-name mapping (31 scripts) is recorded in [`RENAME_MAP.tsv`](RENAME_MAP.tsv). Renames did not alter statistical code or data.

## Structure

```
├── blood/                     ← blood transcriptome analyses
│   ├── GSE98793/  GSE39653/  GSE19738/  GSE251778/   (expression pipelines)
│   ├── GSE61672/  GSE228702/                           (anxiety / SAD×ELA boundary)
│   └── cross_cohort/           ← blood cross-cohort synthesis (CAMERA/FRY, singscore, meta, HC3, cell proxy)
├── brain/                     ← postmortem brain analyses
│   ├── GSE53987/  GSE54563/  GSE80655/
│   └── cross_cohort/           ← brain boundary synthesis
├── methylation/               ← methylomic analyses
│   ├── GSE251780/              (paired RNA–DNAm, promoter v322, cell-adjusted)
│   ├── GSE298549/  GSE88890/  GSE113725/   (external methylation)
├── cross_cohort/              ← cross-transcriptome targeted preranked OXPHOS validation
├── genetics/                  ← FUMA/MAGMA, SMR analyses
├── lookups/                   ← author-deposited result lookups (Erady DGE, Shen EWAS CpGs, DOT1L)
└── figures/                   ← main-text figure reproduction (Figures 1–5)
```

## Datasets and data acquisition

| GSE | Modality / tissue | Analytic N | Manuscript | Raw data source |
|---|---|---|---|---|
| GSE98793 | Blood microarray (GPL570) | 192 | §2.2, Table 1, Fig 1–2 | GEO series matrix |
| GSE39653 | Blood microarray (GPL10558) | 45 | §2.1–2.4, Tables 1–3 | GEO series matrix / deposit |
| GSE19738 | Blood microarray, basal only (GPL6848) | 67 | §2.1–2.4, Tables 1–3 | GEO series matrix |
| GSE251778 | Blood RNA-seq (NovaSeq) | 169 | §2.1–2.6, Fig 1–4, Table 2 | GEO + author-processed count matrix (supplementary of the source paper) |
| GSE61672 | Community blood microarray | 531 / 336 | §2.5, Table S2 | GEO series matrix |
| GSE228702 | Blood 3′-RNA-seq, SAD × ELA | 160 | §2.5, Table S2 | GEO + supplementary |
| GSE53987 | Brain microarray, HIP/BA46/STR | 32/34/32 | §2.2–2.3, Fig 3 | GEO series matrix |
| GSE54563 | Brain microarray, ACC | 50 | §2.2–2.3, Fig 3 | GEO series matrix |
| GSE80655 | Brain RNA-seq, AnCg/nAcc/DLPFC | 48/44/47 | §2.2–2.3, Fig 3 | GEO raw counts |
| GSE251780 | Blood EPIC methylation (paired) | 147 | §2.6, Fig 4–5, Tables S6/S27 | GEO + author-processed matrix (supplementary) |
| GSE298549 | Blood 450K, wLE (external) | 102 | Text S8, Table S28 | GEO (author-processed M-values) |
| GSE88890 | Brain 450K, BA11/BA25 (external) | 75 spec. | Text S8, Table S28 | GEO series matrix + normalised betas |
| GSE113725 | Blood 450K, depression history (external) | 194 | Text S8, Table S28 | GEO series matrix |

All series matrices, count tables, manifests and normalised-beta files are public downloads from the NCBI GEO page of the respective accession (`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE<accession>`); file names used by each script are given in the table below.

## Gene-set definition

GO:0006119 is used in two annotation-derived versions:
- **150-gene (primary, Supplementary Table S7b)**: 137 nuclear + 13 mitochondrial — used by Methods and Table 2.
- **148-gene (frozen v3.22, sensitivity)**: `org.Hs.eg.db 3.22.0` — used by `blood/GSE251778/GSE251778__blood_rnaseq__weighted_camera_sensitivity.R` (non-primary sensitivity reconstruction, not a reported number). Version difference: S7b-only AFG1L/CHCHD10/COA6/NDUFA4/NIPSNAP2/TAFAZZIN; v3.22-only C2ORF69/COXFA4/ETFRF1/MTCH2.

## Inputs / Outputs per pipeline

Inputs are listed as logical paths under `DATA_ROOT` (public GEO/supplementary downloads) or `WORK_ROOT` / `ARCHIVE_ROOT` (author-produced intermediate tables not distributed with this repository). Outputs are written under `OUTPUT_ROOT` (default `./output/`), relative subpaths as shown. Where one pipeline consumes another's output (e.g. #3 ← #2, #25 ← #26, #29 ← #28), run the producer first.

| # | Script | Inputs | Outputs (under `OUTPUT_ROOT`) |
|---|---|---|---|
| 1 | `lookups/MULTI__blood_expression__published_result_lookup.R` | `WORK_ROOT/transcriptomics/blood/{GSE98793,GSE39653}/results_fixed/*.csv` | `lookups/` |
| 2 | `lookups/GSE251780__blood_methylation__shen_cpg_lookup.py` | `WORK_ROOT/analysis_framework/GSE251780/` + Nature Table 2 (network) | `lookups/` |
| 3 | `lookups/GSE251780__blood_methylation__dot1l_influence_audit.R` | outputs of #2 + `WORK_ROOT/analysis_framework/GSE251780/{data,results}/…` | `lookups/` |
| 4 | `blood/GSE19738/GSE19738__blood_microarray__differential_expression.R` | `DATA_ROOT`: GSE19738 series matrix (or GEOquery download) | `blood/GSE19738/{results,unified_tables}` |
| 5 | `blood/GSE228702/GSE228702__blood_rnaseq__sad_ela_boundary.R` | `WORK_ROOT/analysis_framework/GSE228702/data/` | `blood/GSE228702/{results,figures,logs}` |
| 6 | `blood/GSE251778/GSE251778__blood_rnaseq__differential_expression.R` | `WORK_ROOT/transcriptomics/blood/GSE251778/data/GSE251778_Mokhtari_rna_counts.csv.gz` | `blood/GSE251778/results_fixed/` |
| 7 | `blood/GSE251778/GSE251778__blood_rnaseq__weighted_camera_sensitivity.R` | `WORK_ROOT/transcriptomics/blood/GSE251778/` (counts + metadata) | `blood/GSE251778/cross_cohort_out/` |
| 8 | `blood/GSE39653/GSE39653__blood_microarray__differential_expression.R` | `DATA_ROOT`: GSE39653 series matrix / manual RData / cleaned rds (script searches candidates) | `blood/GSE39653/results_fixed/` |
| 9 | `blood/GSE61672/GSE61672__blood_microarray__anxiety_boundary.R` | `WORK_ROOT/analysis_framework/GSE61672/data/GSE61672_series_matrix.txt.gz` | `blood/GSE61672/{data,results,figures,logs}` |
| 10 | `blood/GSE98793/GSE98793__blood_microarray__differential_expression.R` | `DATA_ROOT`: GSE98793 series matrix (GSE98793_analysis/; ~125 MB) | `blood/GSE98793/results_fixed/` |
| 11 | `brain/GSE53987/GSE53987__brain_microarray__differential_expression.R` | `DATA_ROOT`: GSE53987_series_matrix.txt | `brain/GSE53987/results_fixed/` |
| 12 | `brain/GSE54563/GSE54563__brain_microarray__differential_expression.R` | `DATA_ROOT`: GSE54563_series_matrix.txt | `brain/GSE54563/results_fixed/` |
| 13 | `brain/GSE80655/GSE80655__brain_rnaseq__differential_expression.R` | `DATA_ROOT`: GSE80655_raw_counts_GRCh38.p13_NCBI.tsv (+ optional series matrix, metadata) | `brain/GSE80655/results_fixed/` |
| 14 | `brain/cross_cohort/MULTI__brain_oxphos__boundary_sensitivity.R` | `WORK_ROOT/` (MOFA inputs) + `ARCHIVE_ROOT/outputs/Claude_review_2026-09-30/GO0006119_membership.tsv` | `brain/brain_boundary/` |
| 15 | `blood/cross_cohort/MULTI__blood_oxphos__singscore_meta_sensitivity.R` | `WORK_ROOT/` (multi-cohort matrices) + `WORK_ROOT/analysis_framework/results/{predefined_OXPHOS_across_transcriptomes.tsv,blood_OXPHOS_core_gene_cross_dataset.tsv}` | `cross_cohort/results_strict/` |
| 16 | `blood/cross_cohort/MULTI__blood_oxphos__voom_hc3_sensitivity.R` | `WORK_ROOT/analysis_framework/{results_strict,sensitivity_2026-10-04}` + `WORK_ROOT/transcriptomics/blood/{GSE98793,GSE251778,GSE251780}/…` | `cross_cohort/sensitivity_2026-10-04/` |
| 17 | `blood/cross_cohort/MULTI__blood_signatures__cell_proxy_validation.R` | `WORK_ROOT/analysis_framework/sensitivity_2026-10-04/erady_validation/` | `cross_cohort/signatures/` |
| 18 | `blood/cross_cohort/MULTI__oxphos__camera_fry_correlation_audit.R` | `WORK_ROOT/analysis_framework/results_strict/…` + `ARCHIVE_ROOT/outputs/MDD_Anxiety_online_followup_2026-09-23/GWAS_results/01_结构层_FUMA_MAGMA/gene_level_cross_trait_all.tsv.gz` | `cross_cohort/Claude_review_2026-09-30/` |
| 19 | `cross_cohort/MULTI__transcriptome_oxphos__targeted_preranked_tests.R` | `WORK_ROOT/` (cross-cohort expression matrices / metadata) | `cross_cohort/` |
| 20 | `figures/GSE251780__blood_methylation__figure5_promoter_v322.R` | `ARCHIVE_ROOT/outputs/manuscript_v6_12_optionB_2026-10-13/` | `fig_genetics/script1_v322_fig5/` |
| 21 | `figures/MULTI__figures__manuscript_v6_12.R` | `ARCHIVE_ROOT/` (result tables; Figure 1 source from `…/manuscript_v6_12_optionB_2026-10-13/figure_source_data/`) | `fig_genetics/script2_manuscript_v6_12/` |
| 22 | `genetics/MULTI__genetics__fuma_mdd_anxiety.py` | `ARCHIVE_ROOT/work/fuma_compare_2026-09-23/{mdd,anxiety_maf001}` | `fig_genetics/script3_fuma_cross_trait/` |
| 23 | `genetics/MULTI__genetics__smr_twas_result_tables.py` | `ARCHIVE_ROOT/outputs/{MDD_Anxiety_online_followup_2026-09-23,FUMA_MDD_Anxiety_cross_trait_2026-09-23}` | `fig_genetics/script4_smr_twas_integrated/` |
| 24 | `methylation/GSE113725/GSE113725__blood_methylation__depression_history.R` | `DATA_ROOT/GSE113725/GSE113725_series_matrix.txt` (1.5 GB) + `DATA_ROOT/GSE88890/GSE88890_RAW/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz` | `methylation/GSE113725/results_fixed/` |
| 25 | `methylation/GSE251780/GSE251780__blood_methylation__cell_adjusted_sensitivity.R` | `WORK_ROOT/analysis_framework/GSE251780/data/{GSE251780_Mokhtari_MatrixProcessed_AverageAndPvalue.txt.gz,GSE251780_metadata.csv}` | `methylation/GSE251780/{results,figures}` |
| 26 | `methylation/GSE251780/GSE251780__blood_methylation__paired_ewas_rna_dnam.R` | `WORK_ROOT/analysis_framework/GSE251780/data/` + RNA counts (`WORK_ROOT/transcriptomics/blood/GSE251778/…`) | `methylation/GSE251780/{results,figures,logs}` |
| 27 | `methylation/GSE251780/GSE251780__blood_methylation__promoter_planA_historical.R` | `WORK_ROOT/analysis_framework/GSE251780/` (+ MitoCarta source xls under the plan-A output `sources/`) | `methylation/secondary_planA_2026-10-05/` |
| 28 | `methylation/GSE251780/GSE251780__blood_methylation__promoter_v322_recalculation.R` | `ARCHIVE_ROOT/outputs/secondary_planA_2026-10-05/` + `WORK_ROOT/analysis_framework/GSE251780/` | `methylation/recalc_v322/` |
| 29 | `methylation/GSE251780/GSE251780__blood_methylation__promoter_v322_table_assembly.R` | outputs of #28 + `ARCHIVE_ROOT/outputs/manuscript_v6_12_optionB_2026-10-13/github_repository_draft/` | `methylation/github_repo_GSE251780/` |
| 30 | `methylation/GSE298549/GSE298549__blood_methylation__weighted_life_events.R` | `WORK_ROOT/analysis_framework/GSE201016/data/{GSE298549_mval_MDD_cases.txt.gz,GSE298549_metadata.csv}` | `methylation/GSE298549/{results,figures,logs}` |
| 31 | `methylation/GSE88890/GSE88890__brain_methylation__suicide_mdd.R` | `DATA_ROOT/GSE88890/{GSE88890_series_matrix.txt,GSE88890_normalisedBetas.csv,RAW/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz}` | `methylation/GSE88890/results_fixed/` |

## Changelog (2026-10-21)

- **Machine-specific paths replaced by placeholder roots** (`DATA_ROOT` / `WORK_ROOT` / `ARCHIVE_ROOT` / `OUTPUT_ROOT`, defined at the top of each script); no statistical code, threshold, gene set or result value was altered.
- **Output directory standardised** to `OUTPUT_ROOT` (default relative `./output/`, git-ignored); output file names unchanged.
- **Bioconductor compatibility fix**: `methylation/GSE113725/...` annotation loading updated to `minfi::getAnnotation()` (current Bioconductor object type; same pattern already used by the GSE298549 script).
- **Verification**: all 31 scripts were run to completion end-to-end before path normalisation (round-1 run-through, 2026-10-21); a full syntax check of all 31 scripts was re-run after normalisation. No full re-execution was performed after normalisation.

## Known issues

1. **`weighted_camera` limma parameter**: in `blood/GSE251778/GSE251778__blood_rnaseq__weighted_camera_sensitivity.R`, `camera(..., weighted=FALSE)` is silently ignored by the installed limma version, so the "unweighted" numbers are identical to the weighted run. This is a statistical caveat to revisit; the script was not altered.
2. **`boundary_sensitivity` input source**: `brain/cross_cohort/MULTI__brain_oxphos__boundary_sensitivity.R` reads the historical MOFA inputs under `WORK_ROOT` per its original design — it reproduces the archived boundary analysis and does **not** chain from fresh DE runs of the other `brain/` scripts in this repository.
3. **Placeholder paths**: the scripts are archived analysis code; the placeholder roots must be edited and inputs supplied before re-running (see [Quick start](#quick-start)). A few scripts reference an optional developer R-library path as `<DEV_R_LIBRARY_PATH>`; on machines without such a library this is a no-op.

## Environment

R 4.5.x; packages include limma/edgeR (voom, CAMERA, FRY), DESeq2, fgsea, singscore, clusterProfiler, org.Hs.eg.db (3.22.0 for gene-set queries), AnnotationDbi, data.table, ggplot2, sandwich/matrixStats, plus minfi/EpiDISH/missMethyl and the Illumina 450K/EPIC annotation packages (methylation pipelines). Python (genetics/lookups): pandas, numpy, matplotlib, scipy, statsmodels, tabulate.

## License

MIT License — see [`LICENSE`](LICENSE). The scripts and documentation in this repository are shared for method transparency and reproducibility. Participant-level data and raw GEO/GWAS inputs remain the property of their respective depositors; this release does not redistribute them.

Manuscript and supplement documents are maintained separately (not part of this scripts repository).

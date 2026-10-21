# Archive Disposition Record — 2026-10-21

Purpose: record every non-trivial decision applied while assembling this release archive (`/Users/tan/Documents/MDD_OXPHOS_release_2026-10-21`) from the working packages. All changes below are **copy/annotate operations**; no analysis was re-run and no manuscript number was changed.

---

## 1. Table 2 — GSE251778 row restored to the manuscript primary

- **Problem found (2026-10-16):** `results/cross_cohort/Table2_correlation_comparisons_11.tsv` had its GSE251778 row silently switched to the N2 rebuild values (133 genes, fixed P = 0.0046327866, Source = `N2_GSE251778_weighted_camera_rebuild.tsv`), 28 seconds after the N2 output was written. No script or commit records this switch.
- **Manuscript state:** Methods (150-gene definition), Table 2 in `MDD_OXPHOS_manuscript_MASTER_v2.3.docx`, and Figure 3 source data all consistently use the historical 136-gene / fixed P = 0.0065075776 result.
- **Decision (author-approved):** restore the archive row to the historical values (`blood_rnaseq.tsv`): Genes 136, ρ = 0.2660218574, Fixed_P = 0.0065075776, Estimated_P = 0.4924895126, FRY_P = 0.1820750643, Direction Up, Correction_family unchanged, Source = `blood_rnaseq.tsv`.
- **Action taken:** row replaced in `results/cross_cohort/Table2_correlation_comparisons_11.tsv`; `blood_rnaseq.tsv` (historical CAMERA output, from `manuscript_v6_2_correlation_reassessment_2026-10-06/tables/`) copied into `results/cross_cohort/` as the now-referenced source file.
- **N2 files retained** as non-primary sensitivity reconstruction: `N2_GSE251778_weighted_camera_rebuild.tsv`, `N2_GSE251778_sessionInfo.txt`, `scripts/pathway_scores/N2_GSE251778_weighted_camera_rebuild.R`. The N2 run did not reproduce the historical specification (133 ≠ 136 genes; all P values differ) and does not correspond to the frozen 148-gene list file per checklist G1; it is therefore **not** a basis for any reported number.

## 2. OXPHOS gene-set definitions — 150 vs 148

- **150-gene (S7b):** 137 nuclear + 13 mitochondrial; extracted from Supplement Table S7b → `data/gene_sets/GO0006119_S7b_150.tsv` (Gene + GenomeOfOrigin).
- **148-gene (v3.22):** 135 nuclear + 13 mitochondrial; frozen list already present at `data/gene_sets/GO0006119_official_148.tsv` (149 lines = header + 148).
- **Version diff (case-normalised):** `data/gene_sets/GO0006119_version_diff_150_vs_148.tsv`:
  - S7b only (6): AFG1L, CHCHD10, COA6, NDUFA4, NIPSNAP2, TAFAZZIN
  - v3.22 only (4): C2ORF69, COXFA4, ETFRF1, MTCH2
- **Root cause:** two different `org.Hs.eg.db` releases were queried with GOALL GO:0006119 (137 vs 135 nuclear symbols). Tracked as checklist G2. The manuscript retains the 150-gene primary definition; the 148-gene list is the frozen re-analysis target (G1).

## 3. Coverage numbers (read-only recomputation, 2026-10-21)

Using the GSE251778 count matrix with the v3.22 symbol mapping and mean-count>10 filter:
- 150-definition ∩ matrix = **134** (the historical 136 is not reproducible under v3.22 — the 2-gene shortfall is part of the G2 version issue)
- 148-definition ∩ matrix = **133** (exactly reproduces the N2 rebuild NGenes)
These values were used only for auditing; no manuscript number was changed.

## 4. Content added / standardised during assembly

| Item | From | Note |
|---|---|---|
| `results/cross_cohort/blood_rnaseq.tsv` | `manuscript_v6_2_correlation_reassessment_2026-10-06/tables/` | Historical CAMERA output; now the Table 2 Source for GSE251778 |
| `results/figure_source/Table1_demographics_by_analysis_sample.tsv` | `figure_table_revision_v7_9/source_data/` | Main-text Table 1 source |
| `data/gene_sets/GO0006119_S7b_150.tsv` | Supplement Table S7b (table index 7) | 150 rows + header |
| `data/gene_sets/GO0006119_version_diff_150_vs_148.tsv` | computed | 10 symbols |
| `manuscript/` | `github_repository_draft/master/` + `FINAL_REVISION_CHECKLIST_v7_15.md` | Current versions + changelogs + manifests |
| `scripts/` | `github_repository_draft/scripts/` | `.Rhistory` excluded |
| `results/`, `data/`, `docs/`, `environment/` | `github_repository_draft/` | Unchanged copies |

## 5. Exclusions (not copied)

- `.DS_Store`, `.Rhistory` (working artifacts)
- `github_repository_draft` internal audit documents not needed for reproducibility (e.g. `MISMATCH_report.md`, `PIPELINE_AUDIT_EVALUATION.md`, review matrices) — kept in the draft package, not in this release
- Raw GEO/GWAS downloads and participant-level matrices — access notes only (`data/`)
- Old manuscript versions v1.0–v2.1 (except the two immediate backups above) — full history remains in `github_repository_draft/master/`

## 6. Open items carried forward

- **G-series (checklist G1/G2/G5/G7):** freeze the definitive gene-set definition, record the two `org.Hs.eg.db` version root cause, audit per-script library paths, lock the environment (renv/Docker) + hashes.
- **N2a/N2b:** the N2 rebuild remains incomplete as a formal task (no frozen-list usage, no alias rule, no difference record); decide whether to complete it or archive it permanently as sensitivity.
- **Manuscript-side (not touched in this archive):** Discussion freeze awaiting author review; reference ordering; gene-set wording consistency (150 in Methods vs S7b); Table 1 Age "NR" rows; ethics/AI statements; figure caption (Figure 5 coverage column wording removed in v2.3).
- **Archive-side:** regenerate `SHA256SUMS.txt` after any future change; scrub local paths before public push (16/28 scripts per `scripts/README.md`).

## 7. Provenance of this archive

Built 2026-10-21 from:
- `outputs/manuscript_v6_12_optionB_2026-10-13/github_repository_draft/` (primary working package; scripts/results/data/docs/environment/figures/manuscript)
- `outputs/manuscript_v6_2_correlation_reassessment_2026-10-06/tables/blood_rnaseq.tsv`
- `outputs/manuscript_v6_12_optionB_2026-10-13/figure_table_revision_v7_9/source_data/Table1_demographics_by_analysis_sample.tsv`

Read-only cross-check of v2.2/v2.3: `manuscript/MDD_OXPHOS_v2.2_cross_validation_2026-10-21.md`.

---

## 8. Second restructure — scripts-only GitHub archive (2026-10-21, author instruction)

Per author decision ("GitHub 只需要 scripts；按之前的格式只放最终分析脚本并用不同 GSE 分类，每个 GSE 的脚本给一个 README；不需要 data/results，README 写明输入输出"):

- `data/`, `results/`, `manuscript/`, `docs/`, `environment/`, `figure_source_data/` and the full-content `README.md` were **removed from this folder** (the complete package remains untouched in `outputs/manuscript_v6_12_optionB_2026-10-13/github_repository_draft/`).
- Scripts were reorganised by layer + GSE, mirroring the earlier `transcriptomics/{blood,brain,methylation}/<GSE>/` convention:
  - `blood/GSE98793|GSE39653|GSE19738|GSE251778|GSE61672|GSE228702/` + `blood/cross_cohort/`
  - `brain/GSE53987|GSE54563|GSE80655/` + `brain/cross_cohort/`
  - `methylation/GSE251780|GSE298549|GSE88890|GSE113725/`
  - `cross_cohort/` (TARGET_PRERANK), `genetics/`, `lookups/`, `figures/`
- Each GSE/layer directory now carries a `README.md` stating script role, inputs and outputs (13 GSE + 6 layer/function READMEs + top-level README).
- SHA256 manifest regenerated for the scripts-only tree (56 entries).

---

## 9. FUMA legacy job 773201 archived (2026-10-21, author instruction)

Per author instruction ("旧版图和结果打包到一个路径，补充资料只保存最新版；Downloads 下文件也打包，旧版归档，保存新版"):

- **Legacy archive created**: `github_repository_draft/results/archive/FUMA_legacy_archive_2026-10-21/` (39 files, 87 MB) containing:
  - `jobs/FUMA_job773201.zip` + unpacked `FUMA_job773201/` (from `~/Downloads`)
  - `figures/` — six Desktop PDFs for job 773201 (geneManhattan, geneQQplot, lociPlot, magma_exp_gtex, QQplot, snpAnnotPlot)
  - `work_results/ANX_job773201/` — old FUMA result tables (`ple/work/fuma_compare_2026-09-23/anxiety/`)
- **Why superseded**: job 773201 = anxiety FUMA run with **MAF=0** (unfiltered); job 773422 = corrected anxiety run (**MAF=0.01**) which is the version the supplement references (Text S7).
- **Byte-level comparison 773201 vs 773422**: `magma.genes.out`, `magma.gsa.out`, GTEx tissue output, `leadSNPs.txt` and `GenomicRiskLoci.txt` boundaries are **identical**; only SNP-level annotation tables (`snps.txt`/`IndSigSNPs.txt`/`genes.txt`) differ. Therefore **no supplement table or figure required a numeric change**; Supplementary Figure S2 (gene QQ) renders byte-identically between the two jobs.
- **Policy**: Supplement references latest versions only (MDD = job 773196, anxiety = job 773422). New versions remain in `~/Downloads` (FUMA_job773196/, FUMA_job773422/, plus zips); legacy archive is the designated archive of record with its own SHA256SUMS.txt.

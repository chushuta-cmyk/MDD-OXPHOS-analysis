# MDD–OXPHOS Multi-Omics Re-analysis — Analysis Scripts

Analysis pipeline for *"Reassessing oxidative phosphorylation associations in depression using public transcriptomic, methylomic and genetic data"* (retrospective exploratory re-analysis). This repository contains **only the analysis scripts**; raw data and result tables are **not** distributed. Each GSE-classified directory has its own README stating inputs and outputs.

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

## Datasets

| GSE | Modality / tissue | Analytic N | Manuscript |
|---|---|---|---|
| GSE98793 | Blood microarray (GPL570) | 192 | §2.2, Table 1, Fig 1–2 |
| GSE39653 | Blood microarray (GPL10558) | 45 | §2.1–2.4, Tables 1–3 |
| GSE19738 | Blood microarray, basal only (GPL6848) | 67 | §2.1–2.4, Tables 1–3 |
| GSE251778 | Blood RNA-seq (NovaSeq) | 169 | §2.1–2.6, Fig 1–4, Table 2 |
| GSE61672 | Community blood microarray | 531 / 336 | §2.5, Table S2 |
| GSE228702 | Blood 3′-RNA-seq, SAD × ELA | 160 | §2.5, Table S2 |
| GSE53987 | Brain microarray, HIP/BA46/STR | 32/34/32 | §2.2–2.3, Fig 3 |
| GSE54563 | Brain microarray, ACC | 50 | §2.2–2.3, Fig 3 |
| GSE80655 | Brain RNA-seq, AnCg/nAcc/DLPFC | 48/44/47 | §2.2–2.3, Fig 3 |
| GSE251780 | Blood EPIC methylation (paired) | 147 | §2.6, Fig 4–5, Tables S6/S27 |
| GSE298549 | Blood 450K, wLE (external) | 102 | Text S8, Table S28 |
| GSE88890 | Brain 450K, BA11/BA25 (external) | 75 spec. | Text S8, Table S28 |
| GSE113725 | Blood 450K, depression history (external) | 194 | Text S8, Table S28 |

## Inputs / outputs

- **Inputs**: download raw series/count matrices from GEO (accessions above; links in each per-GSE README). Scripts expect a local `DATA_DIR`; **20 of 31 scripts still contain local absolute paths** — replace with your own paths before running (see pre-release checklist below). Input files are documented per-GSE (file names and GEO source) in each README's *Pipeline files, inputs and outputs* section.
- **Outputs**: each script writes DEG/pathway/score tables (also not distributed here). The manuscript tables/figures were assembled from these outputs; mapping is documented per-GSE in each README.

## Gene-set definition

GO:0006119 is used in two annotation-derived versions:
- **150-gene (primary, Supplementary Table S7b)**: 137 nuclear + 13 mitochondrial — used by Methods and Table 2.
- **148-gene (frozen v3.22, sensitivity)**: `org.Hs.eg.db 3.22.0` — used by `blood/GSE251778/GSE251778__blood_rnaseq__weighted_camera_sensitivity.R` (non-primary sensitivity reconstruction, not a reported number). Version difference: S7b-only AFG1L/CHCHD10/COA6/NDUFA4/NIPSNAP2/TAFAZZIN; v3.22-only C2ORF69/COXFA4/ETFRF1/MTCH2.

## Environment

R 4.5.x; packages include limma/edgeR (voom, CAMERA, FRY), DESeq2, fgsea, singscore, clusterProfiler, org.Hs.eg.db (3.22.0 for gene-set queries), AnnotationDbi, data.table, ggplot2, sandwich/matrixStats. See `ARCHIVE_NOTES_2026-10-21.md` for outstanding environment-locking tasks.

## License

MIT License — see [`LICENSE`](LICENSE). The scripts and documentation in this
repository are shared for method transparency and reproducibility. Participant-level
data and raw GEO/GWAS inputs remain the property of their respective depositors;
this release does not redistribute them.

## Integrity

`SHA256SUMS.txt` at repository root records hashes of every file in this archive. Regenerate after any change:
```bash
find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
```

## Pre-release checklist

1. Scrub local absolute paths (`/Users/tan/...`, machine-specific `.libPaths`) → configurable `DATA_DIR` / `RESULTS_DIR` / `OUT`.
2. Credential / personal-identifier scan (no tokens, no participant-level data).
3. License: MIT (see [`LICENSE`](LICENSE)) — confirm it matches your preference before public release.
4. Exclude non-pipeline files (`.Rhistory`, `.DS_Store`) — already removed from this release.
5. Regenerate `SHA256SUMS.txt`; verify before push.
6. Recommended: Zenodo DOI for the final public snapshot.

Manuscript and supplement documents are maintained separately (not part of this scripts repository).

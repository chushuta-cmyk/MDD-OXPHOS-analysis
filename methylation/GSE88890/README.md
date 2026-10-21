# GSE88890 — Postmortem Brain 450K Methylation (Suicide / MDD, BA11 & BA25; Exploratory)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE88890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88890) |
| **Title** | Methylomic profiling of cortex samples from completed suicide cases implicates a role for PSORS1C3 in major depression and suicide |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina HumanMethylation450 BeadChip |
| **Technology** | DNA methylation array (postmortem brain) |
| **Reference** | Murphy TM, Crawford B, Dempster EL, et al. *Transl Psychiatry*. 2017;7:e989. doi:10.1038/tp.2016.249 |

---

## 2. Sample Design

| Region | Cases / Controls |
|--------|-----------------:|
| BA11 (orbitofrontal cortex) | 20 / 20 |
| BA25 (subgenual cingulate) | 17 / 18 |
| **Total specimens** | **75** from **40 donors** (repeated regions per donor) |

- **Phenotype**: Completed suicide with MDD vs non-psychiatric controls.
- **Retained probes**: 347,724 (after beta → M conversion and QC).

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Exploratory external methylation model (Supplementary Text S8); suicide-brain tissue, not the paired blood MDD promoter analysis |
| **Manuscript sections** | Supplementary Text S8; Table S28; Supplementary Data S17 |

---

## 4. Processing in This Re-analysis

From Supplementary Text S8:

- Processed brain beta values transformed to M values; 347,724 retained probes.
- **Pooled design**: `phenotype + region + gender + age` with `duplicateCorrelation` and subject blocking across the 75 specimens / 40 donors; saved consensus correlation = 0.0663.
- **Region-stratified designs**: BA11 and BA25 separately with `phenotype + gender + age`; default eBayes.

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE88890__brain_methylation__suicide_mdd.R` |
| Required inputs | `GSE88890_series_matrix.txt`, `GSE88890_normalisedBetas.csv`, and the GPL13534 manifest `GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz`, currently referenced from a local `GSE88890_RAW/` folder. None is included here. |
| Main outputs written by the script | `GSE88890_DMP_main_blockedSubject.csv`, `GSE88890_DMP_<region>.csv`, `GSE88890_DMP_annotated_adjP0.05.csv`, `GSE88890_chr6_MHC_PSORS1C3_check.csv`, `GSE88890_pvalue_diagnostics.csv`, volcano PDF, report, and session information. |
| Selected outputs included here | results/by_cohort/GSE88890/ *(results are not included in this scripts-only release; see manuscript archive)* contains the annotated significant DMP table, chr6/MHC check, p-value diagnostics, volcano PDF, and report. Full pooled/region-specific DMP tables are not included in this selected archive. |
| Path/run note | Input and output paths are hard-coded to the original workstation. The GPL13534 manifest is shared with the GSE113725 platform, but this script references its local copy under the GSE88890 folder. |

---

## 5. Key OXPHOS Results in v6.12

- Reported only as an exploratory external methylation model; not pooled with paired blood MDD promoter tests.
- All original P and BH q values are preserved as saved; no re-ran DMR or pathway enrichment for OXPHOS in v6.12.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88890>
- Obtain deposited methylation data from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `data.table, limma`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE88890__brain_methylation__suicide_mdd.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Murphy TM, Crawford B, Dempster EL, et al. Methylomic profiling of cortex samples from completed suicide cases implicates a role for PSORS1C3 in major depression and suicide. *Transl Psychiatry*. 2017;7:e989.

# GSE113725 — Whole-Blood 450K Methylation by Self-Reported Depression History (Exploratory)

## 1. GEO Information

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE113725](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113725) |
| **Title** | DNA methylation and inflammation marker profiles associated with a history of depression |
| **Species** | *Homo sapiens* |
| **Platform** | Illumina HumanMethylation450 BeadChip (GPL13534) |
| **Technology** | DNA methylation array (whole blood) |
| **Reference** | Crawford B, Craig Z, Mansell G, et al. *Hum Mol Genet*. 2018;27(16):2840–2850. doi:10.1093/hmg/ddy199 |

---

## 2. Sample Design

| Group | N | Definition |
|-------|--:|------------|
| Groups 1 + 2 | 98 | Self-reported history of doctor-diagnosed depression (± inflammatory disorder) |
| Groups 3 + 4 | 96 | No depression history (± inflammatory disorder) |
| **Total retained** | **194** | |

- **Tissue**: Whole blood.
- **Retained probes**: 350,811.
- **Phenotype**: Self-reported history of depression (not current MDD diagnosis).
- **Saved inflation diagnostic**: λ = 1.216.

---

## 3. Role in This Project

| Aspect | Description |
|--------|-------------|
| **Role** | Exploratory external methylation model (Supplementary Text S8); self-reported depression history, distinct clinical definition from the paired GSE251780 clinical-MDD analysis |
| **Manuscript sections** | Supplementary Text S8; Table S28; Supplementary Data S18 |
| **Not included** | Blood-cell estimates, array/slide adjustment, and regional DMR analysis were not part of this v6.12 run. |

---

## 4. Processing in This Re-analysis

From Supplementary Text S8:

- M-value design: `group + age + gender`; default limma eBayes.
- Group contrast: 1+2 (depression history) vs 3+4 (no history).

### Pipeline files, inputs and outputs

| Item | Files / location |
|---|---|
| Script | `GSE113725__blood_methylation__depression_history.R` |
| Required inputs | `GSE113725_series_matrix.txt` and the GPL13534 manifest `GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz`. The matrix and manifest paths in the current script are absolute workstation paths; neither file is included here. |
| Main outputs written by the script | `GSE113725_DMP_all.csv`, `GSE113725_DMP_significant_adjP0.05.csv`, `GSE113725_pvalue_diagnostics.csv`, `GSE113725_DMP_annotated_adjP0.05.csv` when annotation is available, `GSE113725_LTB4R2_check.csv` when present, volcano PDF, report, and session information. |
| Selected outputs included here | results/by_cohort/GSE113725/ *(results are not included in this scripts-only release; see manuscript archive)* contains the significant DMP table, p-value diagnostics, LTB4R2 check, and report. The full DMP table and source matrix are not included. |
| Path/run note | The script reads the GPL13534 manifest from a local GSE88890 folder; this is a platform manifest dependency, not a GSE113725 expression/methylation sample file. Output is hard-coded to the original workstation. |

---

## 5. Key OXPHOS Results in v6.12

- Reported only as an exploratory external methylation model; not pooled with paired blood MDD promoter tests.
- All original P and BH q values are preserved as saved.

---

## 6. Data Access

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113725>
- Obtain deposited methylation data from GEO. This repository does not redistribute participant-level data.

---

## Usage

### Prerequisites
- R (>= 4.3) with the following Bioconductor/CRAN packages:
  `IlluminaHumanMethylation450kanno.ilmn12.hg19, data.table, limma`

### Run
From this directory, run the pipeline script(s) in order:

```bash
Rscript ./GSE113725__blood_methylation__depression_history.R
```

> Packaging note: input files (GEO series matrices / count tables) are **not** included in this scripts-only release; download them from GEO using the accession in the Dataset Overview. Some scripts reference workstation-specific absolute paths and were not re-run for packaging.

## 7. References

- Crawford B, Craig Z, Mansell G, et al. DNA methylation and inflammation marker profiles associated with a history of depression. *Hum Mol Genet*. 2018;27(16):2840–2850.

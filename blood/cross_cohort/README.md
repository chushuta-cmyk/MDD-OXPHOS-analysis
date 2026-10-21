# Blood Cross-Cohort Synthesis

Cross-cohort OXPHOS analyses over the blood transcriptome datasets.

| Script | Role |
|---|---|
| `MULTI__blood_oxphos__singscore_meta_sensitivity.R` | Fixed GO:0006119 sensitivity + cross-cohort integration (singscore + REML meta; GSE98793/GSE39653/GSE19738 primary; GSE251778 sensitivity-only) |
| `MULTI__blood_oxphos__voom_hc3_sensitivity.R` | Step 3.1 corrections: GSE251778 CAMERA/FRY with full voom weights; GSE98793 nested models & GSE251780 paired models OLS+HC3; sample-overlap check |
| `MULTI__oxphos__camera_fry_correlation_audit.R` | Microarray CAMERA audit (GSE98793/GSE39653/GSE19738) with GO:0006119 membership export |
| `MULTI__blood_signatures__cell_proxy_validation.R` | Cell-proxy / signature validation across blood cohorts |

Inputs: per-cohort expression matrices (see each GSE README). Outputs: Table 2 CAMERA/FRY values, blood meta, cell-marker sensitivity tables.

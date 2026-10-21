#!/usr/bin/env python3
"""Re-audit and integrate current MDD/Anxiety FUMA, SMR, TWAS and MR results."""

from __future__ import annotations

import json
import math
import re
import gzip
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# --- Path configuration: EDIT these placeholder roots for your machine ----------
DATA_ROOT    = "<DATA_ROOT>"
WORK_ROOT    = "<WORK_ROOT>"
ARCHIVE_ROOT = "<ARCHIVE_ROOT>"
OUTPUT_ROOT  = "./output"

ROOT = Path(ARCHIVE_ROOT)
ONLINE = ROOT / "outputs" / "MDD_Anxiety_online_followup_2026-09-23"
FUMA = ROOT / "outputs" / "FUMA_MDD_Anxiety_cross_trait_2026-09-23"
OUT = Path(OUTPUT_ROOT) / "fig_genetics/script4_smr_twas_integrated"
OUT.mkdir(parents=True, exist_ok=True)


def to_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def bh(p: pd.Series) -> np.ndarray:
    x = pd.to_numeric(p, errors="coerce").to_numpy(float)
    out = np.full(len(x), np.nan)
    ok = np.isfinite(x)
    if ok.any():
        out[ok] = multipletests(x[ok], method="fdr_bh")[1]
    return out


def read_smr(path: Path, trait: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    d = pd.read_csv(path, sep="\t", low_memory=False)
    numeric = [
        "p_SMR", "p_SMR_multi", "p_HEIDI", "nsnp_HEIDI", "b_SMR", "se_SMR",
        "p_GWAS", "p_eQTL", "ProbeChr", "Probe_bp", "topSNP_bp",
    ]
    for c in numeric:
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d["trait"] = trait
    d["symbol"] = d["index"].astype(str).str.strip()
    d["gene_id_clean"] = d["gene_id"].astype(str).str.replace(r"\..*$", "", regex=True)
    n = int(d["p_SMR"].notna().sum())
    threshold = 0.05 / n
    d["p_SMR_BH"] = bh(d["p_SMR"])
    d["SMR_Bonf"] = d["p_SMR"] < threshold
    d["HEIDI_gt_0.01"] = d["p_HEIDI"] > 0.01
    d["HEIDI_gt_0.05"] = d["p_HEIDI"] > 0.05
    d["HEIDI_nsnp_ge_10"] = d["nsnp_HEIDI"] >= 10
    d["strict_primary"] = d["SMR_Bonf"] & d["HEIDI_gt_0.01"] & d["HEIDI_nsnp_ge_10"]
    d["strict_sensitivity"] = d["SMR_Bonf"] & d["HEIDI_gt_0.05"] & d["HEIDI_nsnp_ge_10"]

    # One best passing record per gene. Prefer the stricter HEIDI sensitivity definition,
    # then the primary definition, then the minimum SMR P value.
    d["rank_class"] = np.select(
        [d["strict_sensitivity"], d["strict_primary"]], [0, 1], default=2
    )
    ordered = d.sort_values(["symbol", "rank_class", "p_SMR"], kind="mergesort")
    genes = ordered.groupby("symbol", as_index=False).first()
    any_flags = d.groupby("symbol")[["strict_primary", "strict_sensitivity", "SMR_Bonf"]].any()
    genes = genes.drop(columns=["strict_primary", "strict_sensitivity", "SMR_Bonf"]).merge(
        any_flags.reset_index(), on="symbol", how="left"
    )
    genes = genes.sort_values("p_SMR")

    stats = {
        "trait": trait,
        "tests": n,
        "qtl_datasets": int(d["qtl_name"].nunique()),
        "bonf_threshold": threshold,
        "bonf_records": int(d["SMR_Bonf"].sum()),
        "primary_records": int(d["strict_primary"].sum()),
        "sensitivity_records": int(d["strict_sensitivity"].sum()),
        "primary_genes": int(genes["strict_primary"].sum()),
        "sensitivity_genes": int(genes["strict_sensitivity"].sum()),
        "fdr05_records": int((d["p_SMR_BH"] < 0.05).sum()),
    }
    return d, genes, stats


def read_twas(path: Path) -> tuple[pd.DataFrame, dict]:
    d = pd.read_csv(path, sep="\t")
    d = d.rename(columns={"Gene Symbol": "symbol", "P-value": "p_TWAS"})
    d["p_TWAS"] = pd.to_numeric(d["p_TWAS"], errors="coerce")
    d["Z_TWAS"] = pd.to_numeric(
        d["Additional Info"].astype(str).str.extract(r"Z-score:\s*([+-]?[0-9.eE-]+)")[0],
        errors="coerce",
    )
    d["reported_p_lt_5e5"] = d["p_TWAS"] < 5e-5
    # Exact model denominator is absent from the filtered web export. 2.5e-6 is a
    # conservative provisional threshold corresponding to 20,000 tested genes.
    d["provisional_p_lt_2_5e6"] = d["p_TWAS"] < 2.5e-6
    stats = {
        "rows": int(len(d)),
        "unique_genes": int(d["symbol"].nunique()),
        "methods": sorted(d["Method"].dropna().unique().tolist()),
        "min_p": float(d["p_TWAS"].min()),
        "reported_p_lt_5e5": int(d["reported_p_lt_5e5"].sum()),
        "provisional_p_lt_2_5e6": int(d["provisional_p_lt_2_5e6"].sum()),
        "has_fusion_rows": bool(d["Method"].astype(str).str.contains("FUSION", case=False).any()),
    }
    return d.sort_values("p_TWAS"), stats


def read_multimr(path: Path, trait: str) -> tuple[pd.DataFrame, dict]:
    d = pd.read_csv(path, sep="\t")
    d["trait"] = trait
    d["pval"] = pd.to_numeric(d["pval"], errors="coerce")
    d["p_BH_all"] = bh(d["pval"])
    d["p_BH_within_method"] = np.nan
    for method, idx in d.groupby("method").groups.items():
        d.loc[idx, "p_BH_within_method"] = bh(d.loc[idx, "pval"])
    stats = {
        "trait": trait,
        "tests": int(len(d)),
        "genes": int(d["gene_name"].nunique()),
        "nominal_p_lt_0.05": int((d["pval"] < 0.05).sum()),
        "BH_all_lt_0.05": int((d["p_BH_all"] < 0.05).sum()),
        "BH_within_method_lt_0.05": int((d["p_BH_within_method"] < 0.05).sum()),
    }
    return d.sort_values("pval"), stats


mdd_smr, mdd_smr_gene, mdd_smr_stats = read_smr(
    ONLINE / "SMR_MDD_results" / "trait_eSMR.merged.tsv", "MDD"
)
anx_smr, anx_smr_gene, anx_smr_stats = read_smr(
    ONLINE / "SMR_ANX_results" / "trait_eSMR.merged.tsv", "Anxiety"
)
twas, twas_stats = read_twas(ONLINE / "webTWAS_MDD_results" / "final.out")
anx_twas_path = ONLINE / "webTWAS_ANX_results" / "final.out"
if anx_twas_path.exists():
    anx_twas, anx_twas_stats = read_twas(anx_twas_path)
else:
    anx_twas = pd.DataFrame(columns=twas.columns)
    anx_twas_stats = {"available": False, "status": "running_or_not_downloaded"}
mdd_multimr, mdd_multimr_stats = read_multimr(
    ONLINE / "multiMR_results" / "MDD_multiMR.txt", "MDD"
)
anx_multimr, anx_multimr_stats = read_multimr(
    ONLINE / "multiMR_results" / "ANX_multiMR.txt", "Anxiety"
)

fuma = pd.read_csv(FUMA / "gene_level_cross_trait_all.tsv.gz", sep="\t")
for c in ["MDD_Bonf", "ANX_Bonf", "Shared_conjunction_Bonf"]:
    fuma[c] = to_bool(fuma[c])
fuma["symbol"] = fuma["SYMBOL"].astype(str)
fuma_by_symbol = fuma.groupby("symbol", as_index=True).agg(
    MDD_Bonf=("MDD_Bonf", "any"),
    ANX_Bonf=("ANX_Bonf", "any"),
    P_MDD=("P_MDD", "min"),
    P_ANX=("P_ANX", "min"),
)


def smr_gene_subset(d: pd.DataFrame, prefix: str) -> pd.DataFrame:
    cols = [
        "symbol", "gene_id_clean", "qtl_name", "p_SMR", "p_SMR_BH", "b_SMR",
        "se_SMR", "p_HEIDI", "nsnp_HEIDI", "strict_primary",
        "strict_sensitivity", "SMR_Bonf", "GWAS_LOCUS",
    ]
    out = d[[c for c in cols if c in d.columns]].copy()
    return out.rename(columns={c: f"{prefix}_{c}" for c in out.columns if c != "symbol"})


shared_smr = smr_gene_subset(mdd_smr_gene, "MDD").merge(
    smr_gene_subset(anx_smr_gene, "ANX"), on="symbol", how="outer"
)
for c in ["MDD_strict_primary", "MDD_strict_sensitivity", "ANX_strict_primary", "ANX_strict_sensitivity"]:
    shared_smr[c] = shared_smr[c].eq(True)
shared_smr["shared_primary"] = shared_smr["MDD_strict_primary"] & shared_smr["ANX_strict_primary"]
shared_smr["shared_sensitivity"] = (
    shared_smr["MDD_strict_sensitivity"] & shared_smr["ANX_strict_sensitivity"]
)
shared_smr["FUMA_both_Bonf"] = shared_smr["symbol"].map(
    fuma_by_symbol["MDD_Bonf"] & fuma_by_symbol["ANX_Bonf"]
).eq(True)
shared_smr["FUMA_P_MDD"] = shared_smr["symbol"].map(fuma_by_symbol["P_MDD"])
shared_smr["FUMA_P_ANX"] = shared_smr["symbol"].map(fuma_by_symbol["P_ANX"])

twas_symbols = set(twas["symbol"])
twas_strict_symbols = set(twas.loc[twas["provisional_p_lt_2_5e6"], "symbol"])
mdd_smr_primary_symbols = set(mdd_smr_gene.loc[mdd_smr_gene["strict_primary"], "symbol"])
mdd_smr_sens_symbols = set(mdd_smr_gene.loc[mdd_smr_gene["strict_sensitivity"], "symbol"])
anx_smr_primary_symbols = set(anx_smr_gene.loc[anx_smr_gene["strict_primary"], "symbol"])
anx_smr_sens_symbols = set(anx_smr_gene.loc[anx_smr_gene["strict_sensitivity"], "symbol"])
fuma_both_symbols = set(fuma.loc[fuma["MDD_Bonf"] & fuma["ANX_Bonf"], "symbol"])

twas_smr_overlap = twas.merge(
    smr_gene_subset(mdd_smr_gene, "MDD"), on="symbol", how="left"
)
twas_smr_overlap["FUMA_both_Bonf"] = twas_smr_overlap["symbol"].isin(fuma_both_symbols)

# Directional concordance for the exact GTEx Cortex SMR reference, where available.
cortex = mdd_smr[
    mdd_smr["qtl_name"].astype(str).str.match(r"eQTL_GTEx_Brain_Cortex(?:_v10)?$")
].copy()
cortex = cortex.sort_values("p_SMR").drop_duplicates("symbol")
cortex = cortex[["symbol", "qtl_name", "b_SMR", "p_SMR", "p_HEIDI", "nsnp_HEIDI"]]
cortex = cortex.rename(columns={c: f"{c}_Cortex" for c in cortex.columns if c != "symbol"})
twas_smr_overlap = twas_smr_overlap.merge(cortex, on="symbol", how="left")
twas_smr_overlap["same_direction_exact_cortex"] = np.where(
    twas_smr_overlap["b_SMR_Cortex"].notna() & twas_smr_overlap["Z_TWAS"].notna(),
    np.sign(twas_smr_overlap["b_SMR_Cortex"]) == np.sign(twas_smr_overlap["Z_TWAS"]),
    np.nan,
)
twas_smr_overlap["Cortex_SMR_strict_sensitivity"] = (
    (twas_smr_overlap["p_SMR_Cortex"] < mdd_smr_stats["bonf_threshold"])
    & (twas_smr_overlap["p_HEIDI_Cortex"] > 0.05)
    & (twas_smr_overlap["nsnp_HEIDI_Cortex"] >= 10)
)

all_symbols = sorted(fuma_both_symbols | twas_symbols | mdd_smr_primary_symbols | anx_smr_primary_symbols)
evidence = pd.DataFrame({"symbol": all_symbols})
evidence["FUMA_MDD_Bonf"] = evidence["symbol"].map(fuma_by_symbol["MDD_Bonf"]).eq(True)
evidence["FUMA_ANX_Bonf"] = evidence["symbol"].map(fuma_by_symbol["ANX_Bonf"]).eq(True)
evidence["FUMA_both_Bonf"] = evidence["symbol"].isin(fuma_both_symbols)
evidence["MDD_TWAS_p_lt_5e5"] = evidence["symbol"].isin(twas_symbols)
evidence["MDD_TWAS_provisional_strict"] = evidence["symbol"].isin(twas_strict_symbols)
evidence["MDD_SMR_primary"] = evidence["symbol"].isin(mdd_smr_primary_symbols)
evidence["MDD_SMR_sensitivity"] = evidence["symbol"].isin(mdd_smr_sens_symbols)
evidence["ANX_SMR_primary"] = evidence["symbol"].isin(anx_smr_primary_symbols)
evidence["ANX_SMR_sensitivity"] = evidence["symbol"].isin(anx_smr_sens_symbols)
evidence["MDD_ANX_SMR_shared_sensitivity"] = evidence["MDD_SMR_sensitivity"] & evidence["ANX_SMR_sensitivity"]
evidence["P_MDD_MAGMA"] = evidence["symbol"].map(fuma_by_symbol["P_MDD"])
evidence["P_ANX_MAGMA"] = evidence["symbol"].map(fuma_by_symbol["P_ANX"])
evidence["P_MDD_TWAS"] = evidence["symbol"].map(twas.set_index("symbol")["p_TWAS"])
evidence["P_MDD_SMR"] = evidence["symbol"].map(mdd_smr_gene.set_index("symbol")["p_SMR"])
evidence["P_ANX_SMR"] = evidence["symbol"].map(anx_smr_gene.set_index("symbol")["p_SMR"])


def classify(r: pd.Series) -> str:
    if r["MDD_ANX_SMR_shared_sensitivity"] and r["FUMA_both_Bonf"]:
        return "Cross-trait SMR + FUMA"
    if r["MDD_ANX_SMR_shared_sensitivity"]:
        return "Cross-trait SMR"
    if r["FUMA_both_Bonf"] and (r["MDD_SMR_sensitivity"] or r["ANX_SMR_sensitivity"] or r["MDD_TWAS_p_lt_5e5"]):
        return "FUMA shared + one molecular result"
    if r["MDD_SMR_sensitivity"] or r["MDD_TWAS_p_lt_5e5"]:
        return "MDD molecular result only"
    if r["ANX_SMR_sensitivity"]:
        return "Anxiety molecular result only"
    return "FUMA shared only"


evidence["evidence_category"] = evidence.apply(classify, axis=1)
priority_order = {
    "Cross-trait SMR + FUMA": 0,
    "Cross-trait SMR": 1,
    "FUMA shared + one molecular result": 2,
    "MDD molecular result only": 3,
    "Anxiety molecular result only": 4,
    "FUMA shared only": 5,
}
evidence["sort_order"] = evidence["evidence_category"].map(priority_order)
evidence = evidence.sort_values(["sort_order", "P_MDD_MAGMA", "P_MDD_SMR"], na_position="last").drop(columns="sort_order")

# MR result audit.
forward_mr = pd.read_csv(ONLINE / "GWAShug_MR_results" / "MR_result.csv")
reverse_mr = pd.read_csv(ONLINE / "GWAShug_MR_results" / "reverse_MR_result.csv")
forward_mr["direction"] = "MDD -> Anxiety"
reverse_mr["direction"] = "Anxiety -> MDD"
mr = pd.concat([forward_mr, reverse_mr], ignore_index=True)
pleio = pd.read_csv(ONLINE / "GWAShug_MR_results" / "pleiotropy.csv")
hetero = pd.read_csv(ONLINE / "GWAShug_MR_results" / "heterogeneity.csv")


def audit_anx_frequency_orientation(sample_every: int = 200, max_rows: int = 200000) -> tuple[dict, pd.DataFrame]:
    """Compare a systematic ANX sample with independently oriented MDD A1 frequencies."""
    anx_path = ONLINE / "03_SMR" / "ANX_EUR_SMR.ma.gz"
    mdd_path = ONLINE / "03_SMR" / "MDD_no23andMe_EUR_SMR.ma.gz"
    target: dict[str, tuple[str, str, float]] = {}
    with gzip.open(anx_path, "rt") as f:
        header = f.readline().split()
        ix = {x: i for i, x in enumerate(header)}
        for i, line in enumerate(f):
            if i >= max_rows:
                break
            if i % sample_every == 0:
                r = line.split()
                target[r[ix["SNP"]]] = (r[ix["A1"]], r[ix["A2"]], float(r[ix["freq"]]))
    found: dict[str, tuple[str, str, float]] = {}
    with gzip.open(mdd_path, "rt") as f:
        header = f.readline().split()
        ix = {x: i for i, x in enumerate(header)}
        for line in f:
            r = line.split()
            snp = r[ix["SNP"]]
            if snp in target:
                found[snp] = (r[ix["A1"]], r[ix["A2"]], float(r[ix["freq"]]))
    rows = []
    for snp, (a1, a2, af) in target.items():
        if snp not in found:
            continue
        b1, b2, bf = found[snp]
        if (a1, a2) == (b1, b2):
            relation = "same"
            expected = bf
        elif (a1, a2) == (b2, b1):
            relation = "swapped"
            expected = 1 - bf
        else:
            relation = "other"
            expected = np.nan
        rows.append([snp, a1, a2, af, b1, b2, bf, relation, expected, abs(af - expected) if np.isfinite(expected) else np.nan])
    audit = pd.DataFrame(rows, columns=[
        "SNP", "ANX_A1", "ANX_A2", "ANX_A1_freq", "MDD_A1", "MDD_A2", "MDD_A1_freq",
        "allele_relation", "MDD_freq_oriented_to_ANX_A1", "absolute_frequency_difference",
    ])
    usable = audit[audit["allele_relation"].isin(["same", "swapped"])]
    stats = {
        "sampled_ANX_SNPs": len(target),
        "matched_MDD_SNPs": int(len(audit)),
        "same_allele_orientation": int((audit["allele_relation"] == "same").sum()),
        "swapped_allele_orientation": int((audit["allele_relation"] == "swapped").sum()),
        "other_allele_mapping": int((audit["allele_relation"] == "other").sum()),
        "median_abs_frequency_difference": float(usable["absolute_frequency_difference"].median()),
        "provenance_script_available": False,
    }
    return stats, audit


anx_frequency_stats, anx_frequency_audit = audit_anx_frequency_orientation()

summary = {
    "SMR_MDD": mdd_smr_stats,
    "SMR_Anxiety": anx_smr_stats,
    "SMR_shared_primary_genes": sorted(shared_smr.loc[shared_smr["shared_primary"], "symbol"].tolist()),
    "SMR_shared_sensitivity_genes": sorted(shared_smr.loc[shared_smr["shared_sensitivity"], "symbol"].tolist()),
    "TWAS_MDD": twas_stats,
    "TWAS_Anxiety": anx_twas_stats,
    "TWAS_MDD_overlap_SMR_primary": sorted(twas_symbols & mdd_smr_primary_symbols),
    "TWAS_MDD_overlap_SMR_sensitivity": sorted(twas_symbols & mdd_smr_sens_symbols),
    "TWAS_MDD_exact_Cortex_SMR_sensitivity": sorted(
        twas_smr_overlap.loc[twas_smr_overlap["Cortex_SMR_strict_sensitivity"], "symbol"].tolist()
    ),
    "TWAS_MDD_exact_Cortex_direction_concordant_count": int(
        twas_smr_overlap.loc[twas_smr_overlap["Cortex_SMR_strict_sensitivity"], "same_direction_exact_cortex"].sum()
    ),
    "TWAS_MDD_overlap_FUMA_both": sorted(twas_symbols & fuma_both_symbols),
    "multiMR_MDD": mdd_multimr_stats,
    "multiMR_Anxiety": anx_multimr_stats,
    "MR_forward_Egger_intercept_p": float(pleio["pval"].iloc[0]),
    "MR_forward_IVW_Q_p": float(hetero.loc[hetero["method"].str.contains("Inverse"), "Q_pval"].iloc[0]),
    "ANX_SMR_frequency_orientation_audit": anx_frequency_stats,
    "LDSC_rg_result_available": False,
    "Anxiety_TWAS_result_available": bool(not anx_twas.empty),
}

# Outputs.
mdd_smr.to_csv(OUT / "SMR_MDD_all_tests_with_QC.tsv.gz", sep="\t", index=False)
anx_smr.to_csv(OUT / "SMR_ANX_all_tests_with_QC.tsv.gz", sep="\t", index=False)
mdd_smr_gene.to_csv(OUT / "SMR_MDD_gene_summary.tsv", sep="\t", index=False)
anx_smr_gene.to_csv(OUT / "SMR_ANX_gene_summary.tsv", sep="\t", index=False)
shared_smr.to_csv(OUT / "SMR_cross_trait_gene_comparison.tsv", sep="\t", index=False)
twas.to_csv(OUT / "TWAS_MDD_reported_results_QC.tsv", sep="\t", index=False)
twas_smr_overlap.to_csv(OUT / "TWAS_MDD_SMR_overlap.tsv", sep="\t", index=False)
evidence.to_csv(OUT / "current_cross_method_evidence_matrix.tsv", sep="\t", index=False)
mr.to_csv(OUT / "GWAShug_bidirectional_MR_methods.tsv", sep="\t", index=False)
pleio.to_csv(OUT / "GWAShug_forward_MR_pleiotropy.tsv", sep="\t", index=False)
hetero.to_csv(OUT / "GWAShug_forward_MR_heterogeneity.tsv", sep="\t", index=False)
mdd_multimr.to_csv(OUT / "multiMR_MDD_with_FDR.tsv", sep="\t", index=False)
anx_multimr.to_csv(OUT / "multiMR_ANX_with_FDR.tsv", sep="\t", index=False)
anx_frequency_audit.to_csv(OUT / "ANX_SMR_frequency_orientation_audit.tsv", sep="\t", index=False)

if not anx_twas.empty:
    anx_twas.to_csv(OUT / "TWAS_ANX_reported_results_QC.tsv", sep="\t", index=False)
    twas_cross = twas[["symbol", "p_TWAS", "Z_TWAS", "Method"]].rename(
        columns={"p_TWAS": "MDD_p_TWAS", "Z_TWAS": "MDD_Z_TWAS", "Method": "MDD_Method"}
    ).merge(
        anx_twas[["symbol", "p_TWAS", "Z_TWAS", "Method"]].rename(
            columns={"p_TWAS": "ANX_p_TWAS", "Z_TWAS": "ANX_Z_TWAS", "Method": "ANX_Method"}
        ), on="symbol", how="outer"
    )
    twas_cross["both_reported"] = twas_cross["MDD_p_TWAS"].notna() & twas_cross["ANX_p_TWAS"].notna()
    twas_cross["same_direction"] = np.where(
        twas_cross["both_reported"],
        np.sign(twas_cross["MDD_Z_TWAS"]) == np.sign(twas_cross["ANX_Z_TWAS"]),
        np.nan,
    )
    twas_cross.to_csv(OUT / "TWAS_MDD_ANX_cross_trait_comparison.tsv", sep="\t", index=False)
else:
    pd.DataFrame(columns=[
        "symbol", "MDD_p_TWAS", "MDD_Z_TWAS", "MDD_Method",
        "ANX_p_TWAS", "ANX_Z_TWAS", "ANX_Method", "both_reported", "same_direction",
    ]).to_csv(OUT / "TWAS_MDD_ANX_cross_trait_comparison.tsv", sep="\t", index=False)
(OUT / "analysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

# Figure 1: evidence matrix for the most informative genes.
top = evidence[
    evidence["MDD_ANX_SMR_shared_sensitivity"]
    | evidence["MDD_SMR_sensitivity"]
    | evidence["ANX_SMR_sensitivity"]
    | evidence["MDD_TWAS_provisional_strict"]
].copy()
top = top.head(35)
heat_cols = [
    "FUMA_both_Bonf", "MDD_TWAS_provisional_strict", "MDD_SMR_sensitivity",
    "ANX_SMR_sensitivity",
]
labels = ["FUMA both", "MDD TWAS*", "MDD SMR", "ANX SMR"]
mat = top[heat_cols].astype(int).to_numpy()
fig_h = max(5.5, 0.27 * len(top) + 1.7)
fig, ax = plt.subplots(figsize=(7.4, fig_h))
ax.imshow(mat, cmap="Blues", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(len(labels)), labels=labels, rotation=25, ha="right")
ax.set_yticks(range(len(top)), labels=top["symbol"])
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        ax.text(j, i, "Yes" if mat[i, j] else "", ha="center", va="center", fontsize=8,
                color="white" if mat[i, j] else "#333333")
ax.set_title("Current cross-method gene evidence")
ax.set_xlabel("*MDD TWAS provisional P < 2.5e-6; exact model denominator unavailable")
ax.tick_params(axis="both", length=0)
fig.tight_layout()
fig.savefig(OUT / "Figure1_current_cross_method_gene_evidence.png", dpi=220)
fig.savefig(OUT / "Figure1_current_cross_method_gene_evidence.pdf")
plt.close(fig)

# Figure 2: bidirectional MR estimates, explicitly showing instability across methods.
plot_mr = mr[np.isfinite(pd.to_numeric(mr["or"], errors="coerce"))].copy()
plot_mr["label"] = plot_mr["method"].str.replace("Inverse variance weighted ", "IVW ", regex=False)
fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharex=False)
for ax, direction in zip(axes, ["MDD -> Anxiety", "Anxiety -> MDD"]):
    x = plot_mr[plot_mr["direction"] == direction].reset_index(drop=True)
    y = np.arange(len(x))
    ax.errorbar(
        x["or"], y,
        xerr=[x["or"] - x["or_lci95"], x["or_uci95"] - x["or"]],
        fmt="o", color="#2166AC", ecolor="#8DA0CB", capsize=3,
    )
    ax.axvline(1, color="#666666", linestyle="--", linewidth=1)
    ax.set_yticks(y, labels=x["label"])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_title(direction)
    ax.set_xlabel("Odds ratio (log scale)")
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.7)
fig.suptitle("GWAShug bidirectional MR estimates")
fig.tight_layout()
fig.savefig(OUT / "Figure2_GWAShug_bidirectional_MR.png", dpi=220)
fig.savefig(OUT / "Figure2_GWAShug_bidirectional_MR.pdf")
plt.close(fig)

print(json.dumps(summary, ensure_ascii=False, indent=2))

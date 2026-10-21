#!/usr/bin/env python3
"""Cross-trait comparison of FUMA/MAGMA results for MDD and anxiety.

Primary shared-signal test: conjunction P = max(P_MDD, P_ANX), which requires
evidence in both traits. ACAT, Fisher, and sample-size-weighted Stouffer are
reported only as exploratory any-trait/combined statistics because participant
overlap between the source GWAS is not available here.
"""

from __future__ import annotations

import configparser
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues, fisher_exact, hypergeom, norm, spearmanr


ROOT = Path("/Users/tan/Documents/Codex/2026-09-19/ple")
MDD = ROOT / "work/fuma_compare_2026-09-23/mdd"
ANX = ROOT / "work/fuma_compare_2026-09-23/anxiety_maf001"
ANX_OLD = ROOT / "outputs/manuscript_v6_12_optionB_2026-10-13/github_repository_draft/results/archive/FUMA_legacy_archive_2026-10-21/work_results/ANX_job773201"
OUT = ROOT / "outputs/FUMA_MDD_Anxiety_cross_trait_2026-09-23"
OUT.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "mdd": "#3366A6",
    "anx": "#D8872F",
    "shared": "#7A4EAB",
    "neutral": "#B7BDC5",
    "dark": "#263238",
    "oxphos": "#2A9D8F",
}


def read_ws(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", comment="#", engine="python")


def bh(p: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    out = np.empty(n)
    out[order] = adjusted
    return out


def acat_two(p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    p1 = np.clip(np.asarray(p1, dtype=float), 1e-15, 1 - 1e-15)
    p2 = np.clip(np.asarray(p2, dtype=float), 1e-15, 1 - 1e-15)
    t = 0.5 * np.tan((0.5 - p1) * np.pi) + 0.5 * np.tan((0.5 - p2) * np.pi)
    out = 0.5 - np.arctan(t) / np.pi
    tiny = np.maximum(p1, p2) < 1e-12
    out[tiny] = np.minimum(1.0, 2.0 / (1.0 / p1[tiny] + 1.0 / p2[tiny]))
    return np.clip(out, 0.0, 1.0)


def config_values(path: Path) -> dict:
    cp = configparser.ConfigParser()
    cp.read(path)
    return {
        "title": cp.get("jobinfo", "title"),
        "N": cp.getint("params", "N"),
        "MAF": cp.getfloat("params", "MAF"),
        "population": cp.get("params", "pop"),
        "reference": cp.get("params", "refpanel"),
        "leadP": cp.getfloat("params", "leadP"),
        "r2_independent": cp.getfloat("params", "r2"),
        "r2_lead": cp.getfloat("params", "r2_2"),
        "merge_distance_kb": cp.getint("params", "mergeDist"),
        "positional_window_kb": cp.getint("posMap", "posMapWindowSize"),
        "genome_build": "GRCh38" if cp.getint("inputfiles", "snp2genegrch38") else "GRCh37",
    }


def table_or(df: pd.DataFrame, a_col: str, b_col: str) -> tuple[float, float]:
    a = df[a_col].astype(bool)
    b = df[b_col].astype(bool)
    tab = [[int((a & b).sum()), int((a & ~b).sum())],
           [int((~a & b).sum()), int((~a & ~b).sum())]]
    odds, p = fisher_exact(tab, alternative="greater")
    return float(odds), float(p)


def add_gene_statistics(mdd: pd.DataFrame, anx: pd.DataFrame, n_mdd: int, n_anx: int) -> pd.DataFrame:
    keep = ["GENE", "SYMBOL", "CHR", "START", "STOP", "NSNPS", "ZSTAT", "P"]
    g = mdd[keep].merge(anx[keep], on="GENE", how="inner", suffixes=("_MDD", "_ANX"))
    g["SYMBOL"] = g["SYMBOL_MDD"].where(g["SYMBOL_MDD"].notna(), g["SYMBOL_ANX"])
    g["P_MDD_BH"] = bh(g["P_MDD"])
    g["P_ANX_BH"] = bh(g["P_ANX"])
    thr_mdd = 0.05 / len(mdd)
    thr_anx = 0.05 / len(anx)
    thr_common = 0.05 / len(g)
    g["MDD_Bonf"] = g["P_MDD"] < thr_mdd
    g["ANX_Bonf"] = g["P_ANX"] < thr_anx
    g["MDD_nominal"] = g["P_MDD"] < 0.05
    g["ANX_nominal"] = g["P_ANX"] < 0.05
    g["P_conjunction"] = np.maximum(g["P_MDD"], g["P_ANX"])
    g["P_conjunction_Bonf"] = np.minimum(g["P_conjunction"] * len(g), 1.0)
    g["P_ACAT_exploratory"] = acat_two(g["P_MDD"].values, g["P_ANX"].values)
    fisher_p = []
    for p1, p2 in zip(g["P_MDD"], g["P_ANX"]):
        fisher_p.append(combine_pvalues([max(p1, 1e-300), max(p2, 1e-300)], method="fisher").pvalue)
    g["P_Fisher_exploratory"] = fisher_p
    w1, w2 = math.sqrt(n_mdd), math.sqrt(n_anx)
    g["Z_Stouffer_exploratory"] = (w1 * g["ZSTAT_MDD"] + w2 * g["ZSTAT_ANX"]) / math.sqrt(w1*w1 + w2*w2)
    g["P_Stouffer_exploratory"] = norm.sf(g["Z_Stouffer_exploratory"])
    g["Shared_conjunction_Bonf"] = g["P_conjunction"] < thr_common
    g["Evidence_tier"] = np.select(
        [
            g["MDD_Bonf"] & g["ANX_Bonf"],
            g["Shared_conjunction_Bonf"],
            (g["MDD_Bonf"] & g["ANX_nominal"]) | (g["ANX_Bonf"] & g["MDD_nominal"]),
            g["MDD_nominal"] & g["ANX_nominal"],
        ],
        ["Bonferroni in both", "Conjunction Bonferroni", "Bonferroni in one + nominal in other", "Nominal in both"],
        default="No cross-trait support",
    )
    cols = [
        "GENE", "SYMBOL", "CHR_MDD", "START_MDD", "STOP_MDD",
        "NSNPS_MDD", "ZSTAT_MDD", "P_MDD", "P_MDD_BH",
        "NSNPS_ANX", "ZSTAT_ANX", "P_ANX", "P_ANX_BH",
        "MDD_Bonf", "ANX_Bonf", "MDD_nominal", "ANX_nominal",
        "P_conjunction", "P_conjunction_Bonf", "Shared_conjunction_Bonf",
        "P_ACAT_exploratory", "P_Fisher_exploratory",
        "Z_Stouffer_exploratory", "P_Stouffer_exploratory", "Evidence_tier",
    ]
    return g[cols].sort_values(["P_conjunction", "P_ACAT_exploratory"]), thr_mdd, thr_anx, thr_common


def overlap_loci(mdd: pd.DataFrame, anx: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, a in anx.iterrows():
        hits = mdd[(mdd["chr"] == a["chr"]) & (mdd["start"] <= a["end"]) & (mdd["end"] >= a["start"])]
        for _, m in hits.iterrows():
            lo = max(int(m["start"]), int(a["start"]))
            hi = min(int(m["end"]), int(a["end"]))
            rows.append({
                "chr": int(a["chr"]),
                "MDD_locus": int(m["GenomicLocus"]), "MDD_lead_rsID": m["rsID"],
                "MDD_start": int(m["start"]), "MDD_end": int(m["end"]), "MDD_P": m["p"],
                "ANX_locus": int(a["GenomicLocus"]), "ANX_lead_rsID": a["rsID"],
                "ANX_start": int(a["start"]), "ANX_end": int(a["end"]), "ANX_P": a["p"],
                "overlap_start": lo, "overlap_end": hi, "overlap_bp": hi - lo + 1,
                "same_lead_SNP": str(m["rsID"]) == str(a["rsID"]),
            })
    return pd.DataFrame(rows).sort_values(["chr", "overlap_start"]) if rows else pd.DataFrame()


def locus_permutation_p(mdd: pd.DataFrame, anx: pd.DataFrame, observed: int, n_perm: int = 20000) -> tuple[float, float, float, np.ndarray]:
    # GRCh37 autosome lengths; chromosome-preserving random placement.
    chr_len = {1:249250621,2:243199373,3:198022430,4:191154276,5:180915260,6:171115067,
               7:159138663,8:146364022,9:141213431,10:135534747,11:135006516,12:133851895,
               13:115169878,14:107349540,15:102531392,16:90354753,17:81195210,18:78077248,
               19:59128983,20:63025520,21:48129895,22:51304566}
    by_chr = {int(c): x[["start", "end"]].to_numpy(dtype=int) for c, x in mdd.groupby("chr")}
    rng = np.random.default_rng(20260923)
    counts = np.zeros(n_perm, dtype=int)
    for _, row in anx.iterrows():
        c = int(row["chr"])
        if c not in chr_len or c not in by_chr:
            continue
        width = int(row["end"] - row["start"] + 1)
        starts = rng.integers(1, max(2, chr_len[c] - width + 2), size=n_perm)
        ends = starts + width - 1
        intervals = by_chr[c]
        any_hit = np.zeros(n_perm, dtype=bool)
        for ms, me in intervals:
            any_hit |= (starts <= me) & (ends >= ms)
        counts += any_hit
    p = (1 + int((counts >= observed).sum())) / (n_perm + 1)
    return float(p), float(counts.mean()), float(np.quantile(counts, 0.95)), counts


def mapped_genes_for_loci(genes: pd.DataFrame, locus_ids: set[int]) -> str:
    if not locus_ids:
        return ""
    ids = genes["GenomicLocus"].astype(str).str.split(";")
    mask = ids.apply(lambda vals: any(v.isdigit() and int(v) in locus_ids for v in vals))
    return ";".join(sorted(set(genes.loc[mask, "symbol"].dropna().astype(str))))


def genes_by_locus(genes: pd.DataFrame) -> dict[int, str]:
    """Expand FUMA semicolon-delimited locus membership into a locus→genes map."""
    out: dict[int, set[str]] = {}
    for _, row in genes.iterrows():
        symbol = str(row["symbol"])
        for token in str(row["GenomicLocus"]).split(";"):
            token = token.strip()
            if token.isdigit():
                out.setdefault(int(token), set()).add(symbol)
    return {k: ";".join(sorted(v)) for k, v in out.items()}


def pathway_table(mdd: pd.DataFrame, anx: pd.DataFrame) -> tuple[pd.DataFrame, float, float, float]:
    p = mdd.merge(anx, on="FULL_NAME", how="inner", suffixes=("_MDD", "_ANX"))
    thr_mdd, thr_anx, thr_common = 0.05 / len(mdd), 0.05 / len(anx), 0.05 / len(p)
    p["MDD_Bonf"] = p["P_MDD"] < thr_mdd
    p["ANX_Bonf"] = p["P_ANX"] < thr_anx
    p["MDD_FDR"] = bh(p["P_MDD"])
    p["ANX_FDR"] = bh(p["P_ANX"])
    p["P_conjunction"] = np.maximum(p["P_MDD"], p["P_ANX"])
    p["P_conjunction_Bonf"] = np.minimum(p["P_conjunction"] * len(p), 1.0)
    p["P_ACAT_exploratory"] = acat_two(p["P_MDD"].values, p["P_ANX"].values)
    p["Shared_conjunction_Bonf"] = p["P_conjunction"] < thr_common
    return p.sort_values(["P_conjunction", "P_ACAT_exploratory"]), thr_mdd, thr_anx, thr_common


def save_gene_plot(g: pd.DataFrame, path: Path) -> None:
    x = -np.log10(np.clip(g["P_MDD"], 1e-300, 1))
    y = -np.log10(np.clip(g["P_ANX"], 1e-300, 1))
    colors = np.where(g["MDD_Bonf"] & g["ANX_Bonf"], PALETTE["shared"],
             np.where(g["MDD_Bonf"], PALETTE["mdd"], np.where(g["ANX_Bonf"], PALETTE["anx"], PALETTE["neutral"])))
    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    ax.scatter(x, y, c=colors, s=8, alpha=0.55, linewidths=0)
    ax.axvline(-math.log10(0.05 / 18773), color=PALETTE["mdd"], lw=1, ls="--")
    ax.axhline(-math.log10(0.05 / 17461), color=PALETTE["anx"], lw=1, ls="--")
    lab = g[g["Evidence_tier"] != "No cross-trait support"].head(18)
    for _, r in lab.iterrows():
        ax.annotate(r["SYMBOL"], (-math.log10(max(r["P_MDD"], 1e-300)), -math.log10(max(r["P_ANX"], 1e-300))),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("MDD MAGMA gene −log10(P)")
    ax.set_ylabel("Anxiety MAGMA gene −log10(P)")
    ax.set_title("Cross-trait MAGMA gene associations")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def save_pathway_plot(p: pd.DataFrame, path: Path) -> None:
    x = -np.log10(np.clip(p["P_MDD"], 1e-300, 1))
    y = -np.log10(np.clip(p["P_ANX"], 1e-300, 1))
    ox = p["FULL_NAME"].str.contains("OXIDATIVE_PHOSPHORYLATION|ELECTRON_TRANSPORT|RESPIRATORY_CHAIN", case=False, regex=True)
    sig = p["MDD_Bonf"] | p["ANX_Bonf"] | p["Shared_conjunction_Bonf"]
    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    ax.scatter(x[~sig & ~ox], y[~sig & ~ox], c=PALETTE["neutral"], s=7, alpha=0.35, linewidths=0)
    ax.scatter(x[ox], y[ox], c=PALETTE["oxphos"], s=24, alpha=0.9, label="OXPHOS/ETC sets")
    ax.scatter(x[sig], y[sig], c=PALETTE["shared"], s=26, alpha=0.9, label="Bonferroni in ≥1 analysis")
    lab = p[sig | (p["P_conjunction"] < 1e-4)].head(18)
    for _, r in lab.iterrows():
        label = str(r["FULL_NAME"]).replace("GOCC_", "").replace("GOBP_", "")[:38]
        ax.annotate(label, (-math.log10(max(r["P_MDD"], 1e-300)), -math.log10(max(r["P_ANX"], 1e-300))),
                    fontsize=6.5, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("MDD MAGMA gene-set −log10(P)")
    ax.set_ylabel("Anxiety MAGMA gene-set −log10(P)")
    ax.set_title("Cross-trait gene-set convergence")
    ax.legend(frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def save_tissue_plot(t: pd.DataFrame, path: Path) -> None:
    x = -np.log10(np.clip(t["P_MDD"], 1e-300, 1))
    y = -np.log10(np.clip(t["P_ANX"], 1e-300, 1))
    brain = t["FULL_NAME"].str.startswith("Brain_")
    whole_blood = t["FULL_NAME"].eq("Whole_Blood")
    colors = np.where(whole_blood, "#C23B22", np.where(brain, PALETTE["shared"], PALETTE["neutral"]))
    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    ax.scatter(x, y, c=colors, s=np.where(brain | whole_blood, 36, 18), alpha=0.85, linewidths=0)
    for _, r in t[brain | whole_blood | t["MDD_Bonf"] | t["ANX_Bonf"]].iterrows():
        ax.annotate(str(r["FULL_NAME"]).replace("Brain_", ""),
                    (-math.log10(max(r["P_MDD"], 1e-300)), -math.log10(max(r["P_ANX"], 1e-300))),
                    fontsize=6.5, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("MDD GTEx tissue −log10(P)")
    ax.set_ylabel("Anxiety GTEx tissue −log10(P)")
    ax.set_title("GTEx tissue enrichment comparison")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def save_locus_plot(mdd: pd.DataFrame, anx: pd.DataFrame, hits: pd.DataFrame, perm_counts: np.ndarray, path: Path) -> None:
    chr_len = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,
               141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
               81195210,78077248,59128983,63025520,48129895,51304566]
    offsets = np.cumsum([0] + chr_len[:-1])
    overlapped = set(hits["ANX_locus"].astype(int)) if len(hits) else set()
    exact = set(hits.loc[hits["same_lead_SNP"], "ANX_locus"].astype(int)) if len(hits) else set()
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 6.8), gridspec_kw={"height_ratios": [1.15, 1]})
    ax = axes[0]
    for i, length in enumerate(chr_len):
        if i % 2 == 0:
            ax.axvspan(offsets[i] / 1e9, (offsets[i] + length) / 1e9, color="#F4F5F6", zorder=0)
    mx = np.array([offsets[int(c)-1] + (int(s)+int(e))/2 for c,s,e in zip(mdd["chr"],mdd["start"],mdd["end"])]) / 1e9
    ax.scatter(mx, np.repeat(1.0, len(mx)), s=10, color=PALETTE["mdd"], alpha=0.55, linewidths=0, label=f"MDD loci (n={len(mdd)})")
    ax_x = np.array([offsets[int(c)-1] + (int(s)+int(e))/2 for c,s,e in zip(anx["chr"],anx["start"],anx["end"])]) / 1e9
    anx_colors = [PALETTE["shared"] if int(i) in overlapped else PALETTE["anx"] for i in anx["GenomicLocus"]]
    ax.scatter(ax_x, np.repeat(0.0, len(ax_x)), s=34, color=anx_colors, alpha=0.95, linewidths=0.4, edgecolors="white", label=f"Anxiety loci (n={len(anx)})")
    exact_mask = anx["GenomicLocus"].astype(int).isin(exact).to_numpy()
    if exact_mask.any():
        ax.scatter(ax_x[exact_mask], np.repeat(0.0, exact_mask.sum()), s=90, marker="*", color="#111111", label="Exact shared lead SNP")
    ax.set_yticks([0, 1], ["Anxiety", "MDD"])
    ax.set_ylim(-0.45, 1.45)
    ax.set_xlabel("Genome position (GRCh37; chromosomes concatenated)")
    ax.set_title("Genome-wide risk loci and physical interval overlap")
    centers = (np.array(offsets) + np.array(chr_len)/2) / 1e9
    ax.set_xticks(centers, [str(i) for i in range(1, 23)], fontsize=7)
    ax.set_xlabel("Chromosome")
    ax.legend(frameon=False, loc="upper center", ncol=3, fontsize=8)
    ax.spines[["top", "right", "left"]].set_visible(False)

    ax2 = axes[1]
    bins = np.arange(perm_counts.max() + 2) - 0.5
    ax2.hist(perm_counts, bins=bins, color="#B7BDC5", edgecolor="white")
    observed = len(overlapped)
    ax2.axvline(observed, color=PALETTE["shared"], lw=2.5, label=f"Observed = {observed}")
    ax2.set_xlabel("Anxiety loci overlapping ≥1 MDD locus under chromosome-preserving null")
    ax2.set_ylabel("Permutations")
    ax2.set_title("Exploratory interval-overlap permutation (20,000 replicates)")
    ax2.legend(frameon=False)
    ax2.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    mdd_cfg, anx_cfg, old_cfg = config_values(MDD / "params.config"), config_values(ANX / "params.config"), config_values(ANX_OLD / "params.config")
    mdd_g, anx_g = read_ws(MDD / "magma.genes.out"), read_ws(ANX / "magma.genes.out")
    genes, thr_mdd, thr_anx, thr_common = add_gene_statistics(mdd_g, anx_g, mdd_cfg["N"], anx_cfg["N"])

    a = int(genes["MDD_Bonf"].sum()); b = int(genes["ANX_Bonf"].sum()); k = int((genes["MDD_Bonf"] & genes["ANX_Bonf"]).sum()); M = len(genes)
    gene_overlap_p = float(hypergeom.sf(k - 1, M, a, b))
    gene_or, gene_fisher = table_or(genes, "MDD_Bonf", "ANX_Bonf")
    rho_gene, rho_gene_p = spearmanr(genes["ZSTAT_MDD"], genes["ZSTAT_ANX"])

    mdd_loci = pd.read_csv(MDD / "GenomicRiskLoci.txt", sep="\t")
    anx_loci = pd.read_csv(ANX / "GenomicRiskLoci.txt", sep="\t")
    locus_hits = overlap_loci(mdd_loci, anx_loci)
    anx_overlap_n = int(locus_hits["ANX_locus"].nunique()) if len(locus_hits) else 0
    mdd_overlap_n = int(locus_hits["MDD_locus"].nunique()) if len(locus_hits) else 0
    locus_perm_p, locus_perm_mean, locus_perm_p95, locus_perm_counts = locus_permutation_p(mdd_loci, anx_loci, anx_overlap_n)

    mdd_map = pd.read_csv(MDD / "genes.txt", sep="\t")
    anx_map = pd.read_csv(ANX / "genes.txt", sep="\t")
    map_shared = sorted(set(mdd_map["symbol"].dropna()) & set(anx_map["symbol"].dropna()))
    if len(locus_hits):
        mdd_by = genes_by_locus(mdd_map)
        anx_by = genes_by_locus(anx_map)
        locus_hits["MDD_positional_genes"] = locus_hits["MDD_locus"].map(mdd_by).fillna("")
        locus_hits["ANX_positional_genes"] = locus_hits["ANX_locus"].map(anx_by).fillna("")
        locus_hits["shared_positional_genes"] = locus_hits.apply(
            lambda r: ";".join(sorted(set(str(r["MDD_positional_genes"]).split(";")) & set(str(r["ANX_positional_genes"]).split(";")) - {""})), axis=1)

    mdd_gsa, anx_gsa = read_ws(MDD / "magma.gsa.out"), read_ws(ANX / "magma.gsa.out")
    pathways, pthr_mdd, pthr_anx, pthr_common = pathway_table(mdd_gsa, anx_gsa)
    rho_path, rho_path_p = spearmanr(-np.log10(pathways["P_MDD"]), -np.log10(pathways["P_ANX"]))
    ox_regex = "OXIDATIVE_PHOSPHORYLATION|ELECTRON_TRANSPORT|RESPIRATORY_CHAIN|RESPIRATORY_ELECTRON|ATP_SYNTHESIS_COUPLED|MITOCHONDRIAL_RESPIR"
    ox = pathways[pathways["FULL_NAME"].str.contains(ox_regex, case=False, regex=True)].copy()
    mito = pathways[pathways["FULL_NAME"].str.contains("MITOCHON", case=False, regex=True)].copy()

    mdd_t = read_ws(MDD / "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out")
    anx_t = read_ws(ANX / "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out")
    tissues = mdd_t.merge(anx_t, on="FULL_NAME", suffixes=("_MDD", "_ANX"))
    tthr_mdd, tthr_anx = 0.05 / len(mdd_t), 0.05 / len(anx_t)
    tissues["MDD_Bonf"] = tissues["P_MDD"] < tthr_mdd
    tissues["ANX_Bonf"] = tissues["P_ANX"] < tthr_anx
    tissues["P_conjunction"] = np.maximum(tissues["P_MDD"], tissues["P_ANX"])
    tissues["P_conjunction_Bonf"] = np.minimum(tissues["P_conjunction"] * len(tissues), 1.0)
    tissues["P_ACAT_exploratory"] = acat_two(tissues["P_MDD"].values, tissues["P_ANX"].values)
    rho_tissue, rho_tissue_p = spearmanr(-np.log10(tissues["P_MDD"]), -np.log10(tissues["P_ANX"]))
    tissues = tissues.sort_values("P_conjunction")

    # MAF sensitivity audit: MAGMA result identity plus positional changes.
    old_loci = pd.read_csv(ANX_OLD / "GenomicRiskLoci.txt", sep="\t")
    old_snps = pd.read_csv(ANX_OLD / "snps.txt", sep="\t")
    new_snps = pd.read_csv(ANX / "snps.txt", sep="\t")
    locus_sens = old_loci.merge(anx_loci, on="rsID", suffixes=("_MAF0", "_MAF001"), how="outer", indicator=True)
    locus_sens["start_change_bp"] = locus_sens.get("start_MAF001") - locus_sens.get("start_MAF0")
    locus_sens["end_change_bp"] = locus_sens.get("end_MAF001") - locus_sens.get("end_MAF0")
    locus_sens["nSNP_change"] = locus_sens.get("nSNPs_MAF001") - locus_sens.get("nSNPs_MAF0")
    gene_identical = (ANX_OLD / "magma.genes.out").read_bytes() == (ANX / "magma.genes.out").read_bytes()
    gsa_identical = (ANX_OLD / "magma.gsa.out").read_bytes() == (ANX / "magma.gsa.out").read_bytes()
    tissue_identical = (ANX_OLD / "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out").read_bytes() == (ANX / "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out").read_bytes()
    mdd_leads = set(mdd_loci["rsID"].dropna().astype(str))
    anx_leads = set(anx_loci["rsID"].dropna().astype(str))
    mdd_ind = set(pd.read_csv(MDD / "IndSigSNPs.txt", sep="\t")["rsID"].dropna().astype(str))
    anx_ind = set(pd.read_csv(ANX / "IndSigSNPs.txt", sep="\t")["rsID"].dropna().astype(str))

    # Export tables.
    genes.to_csv(OUT / "gene_level_cross_trait_all.tsv.gz", sep="\t", index=False, compression="gzip")
    genes[genes["Evidence_tier"] != "No cross-trait support"].to_csv(OUT / "prioritized_cross_trait_genes.tsv", sep="\t", index=False)
    genes.head(100).to_csv(OUT / "top100_conjunction_genes.tsv", sep="\t", index=False)
    locus_hits.to_csv(OUT / "overlapping_genomic_risk_loci.tsv", sep="\t", index=False)
    pd.DataFrame({"shared_positional_gene": map_shared}).to_csv(OUT / "shared_FUMA_positional_genes.tsv", sep="\t", index=False)
    pathways.to_csv(OUT / "gene_set_cross_trait_all.tsv.gz", sep="\t", index=False, compression="gzip")
    pathways.head(100).to_csv(OUT / "top100_conjunction_gene_sets.tsv", sep="\t", index=False)
    ox.to_csv(OUT / "OXPHOS_ETC_targeted_gene_sets.tsv", sep="\t", index=False)
    mito.to_csv(OUT / "mitochondrial_gene_sets.tsv", sep="\t", index=False)
    tissues.to_csv(OUT / "GTEx_tissue_cross_trait.tsv", sep="\t", index=False)
    locus_sens.to_csv(OUT / "anxiety_MAF_sensitivity_loci.tsv", sep="\t", index=False)

    summary = {
        "input_parameters": {"MDD": mdd_cfg, "Anxiety_MAF001": anx_cfg, "Anxiety_original": old_cfg},
        "gene_analysis": {
            "MDD_tested": int(len(mdd_g)), "Anxiety_tested": int(len(anx_g)), "common_tested": int(M),
            "MDD_Bonferroni_threshold": thr_mdd, "Anxiety_Bonferroni_threshold": thr_anx,
            "common_conjunction_Bonferroni_threshold": thr_common,
            "MDD_significant_in_common": a, "Anxiety_significant_in_common": b, "significant_in_both": k,
            "hypergeometric_overlap_P": gene_overlap_p, "Fisher_overlap_odds_ratio": gene_or,
            "Fisher_overlap_P": gene_fisher, "Spearman_Z_rho": float(rho_gene), "Spearman_Z_P": float(rho_gene_p),
            "conjunction_Bonferroni_genes": int(genes["Shared_conjunction_Bonf"].sum()),
            "at_least_one_Bonf_other_nominal_including_both": int((((genes["MDD_Bonf"] & genes["ANX_nominal"]) | (genes["ANX_Bonf"] & genes["MDD_nominal"]))).sum()),
            "additional_one_Bonf_other_nominal_excluding_both": int(((((genes["MDD_Bonf"] & genes["ANX_nominal"]) | (genes["ANX_Bonf"] & genes["MDD_nominal"]))) & ~(genes["MDD_Bonf"] & genes["ANX_Bonf"])).sum()),
            "nominal_in_both": int((genes["MDD_nominal"] & genes["ANX_nominal"]).sum()),
        },
        "locus_analysis": {
            "MDD_loci": int(len(mdd_loci)), "Anxiety_loci": int(len(anx_loci)),
            "overlap_pairs": int(len(locus_hits)), "MDD_loci_with_overlap": mdd_overlap_n,
            "Anxiety_loci_with_overlap": anx_overlap_n, "chromosome_preserving_permutation_P": locus_perm_p,
            "permutation_mean_Anxiety_loci_overlapping": locus_perm_mean,
            "permutation_95th_percentile": locus_perm_p95,
            "MDD_positional_genes": int(mdd_map["symbol"].nunique()),
            "Anxiety_positional_genes": int(anx_map["symbol"].nunique()),
            "shared_positional_genes": int(len(map_shared)),
            "shared_lead_SNPs": sorted(mdd_leads & anx_leads),
            "n_shared_lead_SNPs": int(len(mdd_leads & anx_leads)),
            "n_shared_independent_significant_SNPs": int(len(mdd_ind & anx_ind)),
        },
        "gene_set_analysis": {
            "MDD_tested": int(len(mdd_gsa)), "Anxiety_tested": int(len(anx_gsa)), "common_tested": int(len(pathways)),
            "MDD_Bonferroni_significant": int((mdd_gsa["P"] < pthr_mdd).sum()),
            "Anxiety_Bonferroni_significant": int((anx_gsa["P"] < pthr_anx).sum()),
            "shared_conjunction_Bonferroni": int(pathways["Shared_conjunction_Bonf"].sum()),
            "Spearman_neglogP_rho": float(rho_path), "Spearman_neglogP_P": float(rho_path_p),
            "OXPHOS_ETC_sets_examined": int(len(ox)),
            "OXPHOS_ETC_min_MDD_P": float(ox["P_MDD"].min()) if len(ox) else None,
            "OXPHOS_ETC_min_Anxiety_P": float(ox["P_ANX"].min()) if len(ox) else None,
            "OXPHOS_ETC_min_conjunction_P": float(ox["P_conjunction"].min()) if len(ox) else None,
        },
        "tissue_analysis": {
            "tissues_tested": int(len(tissues)), "MDD_Bonferroni_significant": int(tissues["MDD_Bonf"].sum()),
            "Anxiety_Bonferroni_significant": int(tissues["ANX_Bonf"].sum()),
            "conjunction_Bonferroni": int((tissues["P_conjunction_Bonf"] < 0.05).sum()),
            "Spearman_neglogP_rho": float(rho_tissue), "Spearman_neglogP_P": float(rho_tissue_p),
            "Whole_Blood_MDD_P": float(tissues.loc[tissues["FULL_NAME"] == "Whole_Blood", "P_MDD"].iloc[0]),
            "Whole_Blood_Anxiety_P": float(tissues.loc[tissues["FULL_NAME"] == "Whole_Blood", "P_ANX"].iloc[0]),
        },
        "MAF_sensitivity": {
            "old_MAF": old_cfg["MAF"], "new_MAF": anx_cfg["MAF"],
            "MAGMA_gene_file_identical": gene_identical, "MAGMA_gene_set_file_identical": gsa_identical,
            "GTEx_tissue_file_identical": tissue_identical,
            "old_candidate_snps": int(len(old_snps)), "new_candidate_snps": int(len(new_snps)),
            "candidate_snps_removed": int(len(set(old_snps["uniqID"]) - set(new_snps["uniqID"]))),
            "lead_SNP_sets_identical": set(old_loci["rsID"]) == set(anx_loci["rsID"]),
            "loci_with_boundary_or_SNP_count_change": int(((locus_sens["start_change_bp"].fillna(0) != 0) | (locus_sens["end_change_bp"].fillna(0) != 0) | (locus_sens["nSNP_change"].fillna(0) != 0)).sum()),
        },
        "interpretation_guardrails": [
            "MDD and anxiety are different phenotypes; overlap is cross-trait convergence, not conventional replication.",
            "Unknown participant overlap can inflate correlations and unadjusted Stouffer/Fisher evidence.",
            "Conjunction max-P is the primary conservative shared-association statistic.",
            "MAGMA Z statistics measure gene association evidence, not biological effect direction.",
            "Formal genetic sharing requires bivariate LDSC; causal prioritization requires TWAS/SMR/colocalization.",
        ],
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    rows = []
    for section, values in summary.items():
        if isinstance(values, dict):
            for metric, value in values.items():
                if not isinstance(value, dict):
                    rows.append({"section": section, "metric": metric, "value": value})
    pd.DataFrame(rows).to_csv(OUT / "analysis_summary.tsv", sep="\t", index=False)

    save_gene_plot(genes, OUT / "Figure1_MAGMA_gene_cross_trait.png")
    save_pathway_plot(pathways, OUT / "Figure2_MAGMA_gene_set_cross_trait.png")
    save_tissue_plot(tissues, OUT / "Figure3_GTEx_tissue_cross_trait.png")
    save_locus_plot(mdd_loci, anx_loci, locus_hits, locus_perm_counts, OUT / "Figure4_genomic_locus_overlap.png")

    # Concise reproducible text report.
    top_shared = genes[genes["Evidence_tier"] != "No cross-trait support"].head(20)[["SYMBOL", "P_MDD", "P_ANX", "P_conjunction", "Evidence_tier"]]
    top_sets = pathways.head(15)[["FULL_NAME", "P_MDD", "P_ANX", "P_conjunction"]]
    report = f"""# MDD–Anxiety FUMA/MAGMA cross-trait analysis

## Scope and primary analysis

This analysis compares MDD job 773196 with the corrected Anxiety MAF=0.01 job 773422. It is a cross-trait convergence analysis rather than a conventional meta-analysis of identical phenotypes. The primary shared-association statistic is the conservative conjunction P value, `max(P_MDD, P_Anxiety)`. ACAT, Fisher, and sample-size-weighted Stouffer results are exploratory because participant overlap is unknown.

## Parameter and sensitivity validation

- MDD: N={mdd_cfg['N']:,}, MAF={mdd_cfg['MAF']}, {mdd_cfg['genome_build']}, {mdd_cfg['population']}.
- Anxiety: N={anx_cfg['N']:,}, MAF={anx_cfg['MAF']}, {anx_cfg['genome_build']}, {anx_cfg['population']}.
- Correcting Anxiety MAF from 0 to 0.01 left the MAGMA gene, gene-set, and GTEx tissue files exactly unchanged.
- Candidate SNPs changed from {len(old_snps):,} to {len(new_snps):,}; {len(set(old_snps['uniqID']) - set(new_snps['uniqID'])):,} candidate SNPs were removed. Lead-SNP sets were unchanged.

## Gene-level convergence

- Common tested genes: {M:,}.
- Bonferroni-significant genes among common genes: MDD {a:,}; Anxiety {b:,}; both {k:,}.
- Significant-gene overlap: hypergeometric P={gene_overlap_p:.3g}; Fisher OR={gene_or:.3g}, P={gene_fisher:.3g}.
- Gene-level MAGMA Z-statistic rank correlation: rho={rho_gene:.3f}, P={rho_gene_p:.3g}. This is descriptive and may partly reflect sample overlap.
- Genes passing conjunction Bonferroni: {int(genes['Shared_conjunction_Bonf'].sum()):,}.
- Genes Bonferroni-significant in both: {k:,}; an additional {int(((((genes['MDD_Bonf'] & genes['ANX_nominal']) | (genes['ANX_Bonf'] & genes['MDD_nominal']))) & ~(genes['MDD_Bonf'] & genes['ANX_Bonf'])).sum()):,} were Bonferroni-significant in one trait and nominal in the other.

Top cross-trait genes:

{top_shared.to_markdown(index=False)}

## Genomic-risk-locus overlap

- MDD loci: {len(mdd_loci)}; Anxiety loci: {len(anx_loci)}.
- Anxiety loci overlapping at least one MDD interval: {anx_overlap_n}; MDD loci overlapping Anxiety: {mdd_overlap_n}; interval pairs: {len(locus_hits)}.
- Chromosome-preserving interval permutation: P={locus_perm_p:.3g}; null mean={locus_perm_mean:.2f}, 95th percentile={locus_perm_p95:.1f} Anxiety loci overlapping.
- Shared positional genes: {len(map_shared)}. Physical interval overlap alone does not establish a shared causal variant; LD-aware colocalization is required.
- Exact shared lead SNPs: {', '.join(sorted(mdd_leads & anx_leads)) if (mdd_leads & anx_leads) else 'none'}; shared independent significant SNPs: {len(mdd_ind & anx_ind)}.

## Gene-set convergence

- Gene sets tested: MDD {len(mdd_gsa):,}; Anxiety {len(anx_gsa):,}; common {len(pathways):,}.
- Bonferroni-significant: MDD {(mdd_gsa['P'] < pthr_mdd).sum()}; Anxiety {(anx_gsa['P'] < pthr_anx).sum()}; conjunction {pathways['Shared_conjunction_Bonf'].sum()}.
- Cross-trait gene-set rank correlation: rho={rho_path:.3f}, P={rho_path_p:.3g}.

Top conjunction gene sets:

{top_sets.to_markdown(index=False)}

## OXPHOS and tissue interpretation

- OXPHOS/ETC targeted sets examined: {len(ox)}. Minimum P values: MDD {ox['P_MDD'].min() if len(ox) else float('nan'):.3g}; Anxiety {ox['P_ANX'].min() if len(ox) else float('nan'):.3g}; conjunction {ox['P_conjunction'].min() if len(ox) else float('nan'):.3g}.
- GTEx tissue correlation: rho={rho_tissue:.3f}, P={rho_tissue_p:.3g}.
- Whole blood: MDD P={summary['tissue_analysis']['Whole_Blood_MDD_P']:.3g}; Anxiety P={summary['tissue_analysis']['Whole_Blood_Anxiety_P']:.3g}.

## Interpretation and next validation

The FUMA/MAGMA comparison can establish cross-trait convergence at genes, pathways, tissues, and broad risk intervals. It cannot by itself establish that overlapping intervals contain the same causal variant, that blood OXPHOS is genetically causal, or that Anxiety independently replicates MDD. The next formal steps are bivariate LDSC for genetic correlation and sample-overlap diagnostics, followed by TWAS/SMR-HEIDI and locus-level colocalization in brain and blood eQTL resources. Blood OXPHOS should remain a peripheral/state-associated layer unless these causal tests support mediation.
"""
    (OUT / "RESULTS_AND_INTERPRETATION.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()

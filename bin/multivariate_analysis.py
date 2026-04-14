#!/usr/bin/env python3
"""
multivariate_analysis.py
━━━━━━━━━━━━━━━━━━━━━━━━
CDC 1.0.0 — Multivariate Statistical Analysis Module
Author : Nadeem Khan, INRS-Centre Armand-Frappier Santé-Biotechnologie
GitHub : github.com/nkhan119

Analyses performed on GWAS summary statistics from Cohort_A (EUR+EAS) and
Cohort_B (AFR+AMR+SAS):

  Phase 1  · Multivariate Hypothesis Testing
             MANOVA, Hotelling T², Pillai/Wilks/Roy/Lawley-Hotelling traces,
             permutation-based significance, cross-trait correlation matrix

  Phase 2  · Regression Analysis
             OLS, Ridge, Lasso (alpha via CV), Elastic Net, PCR, PLS,
             effect-size comparison, conditional R², VIF inflation diagnostics

  Phase 3  · Multiple-Testing Correction
             Bonferroni, Holm-Bonferroni, BH-FDR, BY-FDR, q-value estimation
             (Storey), genomic inflation λ_GC, LDSC intercept proxy

  Phase 4  · Causal Inference
             DoWhy/EconML IV estimation (2SLS), propensity-score matching,
             doubly-robust AIPW, sensitivity (Rosenbaum bounds), effect
             heterogeneity forest

Outputs
-------
  results/multivariate/
    ├── phase1_hypothesis/        TSV tables
    ├── phase2_regression/        TSV tables
    ├── phase3_mtc/               TSV tables
    ├── phase4_causal/            TSV tables
    └── Multivariate_Report.html  single-file interactive report
"""

import os
import sys
import gzip
import glob
import json
import argparse
import warnings
import traceback
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import inv
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import (Ridge, Lasso, ElasticNet,
                                   RidgeCV, LassoCV, ElasticNetCV)
from sklearn.decomposition import PCA as skPCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_score, KFold
from sklearn.metrics import r2_score
from sklearn.neighbors import NearestNeighbors
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Companion_modules ──────────────────────────────────
import sys, importlib
_BIN = str(Path(__file__).parent)
if _BIN not in sys.path:
    sys.path.insert(0, _BIN)

try:
    import static_figures as SF
    _HAS_SF = True
except ImportError as _e:
    print(f"  ⚠ static_figures not importable ({_e}) — skipping static plots")
    _HAS_SF = False

try:
    import excel_tables as ET
    _HAS_ET = True
except ImportError as _e:
    print(f"  ⚠ excel_tables not importable ({_e}) — skipping Excel output")
    _HAS_ET = False

warnings.filterwarnings("ignore")

# ── Plotly_theme ────────────────────────────────────────────────────────────
PALETTE = {
    "primary":   "#2C5F8A",
    "secondary": "#E05A3A",
    "tertiary":  "#3AAE7A",
    "quaternary":"#9B59B6",
    "neutral":   "#7F8C8D",
    "bg":        "#FFFFFF",
    "grid":      "rgba(0,0,0,0)",
    "text":      "#1A1A2E",
    "cohort_A":  "#2C5F8A",
    "cohort_B":  "#E05A3A",
}

COHORT_COLORS = {"Cohort_A": PALETTE["cohort_A"], "Cohort_B": PALETTE["cohort_B"]}

PLOTLY_BASE = dict(
    paper_bgcolor=PALETTE["bg"],
    plot_bgcolor=PALETTE["bg"],
    font=dict(family="Inter, Arial, sans-serif", size=13, color=PALETTE["text"]),
    xaxis=dict(showgrid=False, zeroline=False, linecolor="#CCCCCC",
               showline=True, ticks="outside", ticklen=4),
    yaxis=dict(showgrid=False, zeroline=False, linecolor="#CCCCCC",
               showline=True, ticks="outside", ticklen=4),
    margin=dict(l=70, r=40, t=70, b=60),
    legend=dict(bgcolor="rgba(255,255,255,0.9)", borderwidth=1,
                bordercolor="#DDDDDD"),
)

def apply_base(fig, **kwargs):
    """Apply base layout to any Plotly figure."""
    layout = {**PLOTLY_BASE, **kwargs}
    fig.update_layout(**layout)
    return fig

# ── I/O helpers ─────────────────────────────────────────────────────────────

# ~19M rows × 14 cols instead of × 30+ cols cuts peak memory by ~50%.
_KEEP_COLS = {"SNP","CHR","BP","A1","A2","BETA","SE","P","Z",
              "A1FREQ","INFO","N","GENPOS"}

def load_sumstats(pattern: str) -> pd.DataFrame:
    """
    Load all harmonised *_sumstats.tsv.gz files matching a glob pattern.
    Only retains columns that exist in _KEEP_COLS to limit peak RAM.
    Trait name is extracted from the filename after stripping the cohort prefix,
    e.g.  Cohort_A_ldl_cholesterol_sumstats.tsv.gz  →  trait = 'ldl_cholesterol'
    """
    files = sorted(glob.glob(pattern, recursive=True))
    if not files:
        raise FileNotFoundError(f"No sumstat files found: {pattern}")

    frames = []
    for fp in files:
        fname   = Path(fp).name                            # Cohort_A_bmi_sumstats.tsv.gz
        cohort  = "Cohort_A" if "Cohort_A" in fname else "Cohort_B"
        # Strip cohort prefix and _sumstats suffix
        # Cohort_A_ldl_cholesterol_sumstats.tsv.gz → ldl_cholesterol
        stem    = fname.replace(".gz","").replace(".tsv","")  # Cohort_A_bmi_sumstats
        stem    = stem.replace(f"{cohort}_","",1)             # bmi_sumstats
        trait   = stem.replace("_sumstats","")                # bmi

        with gzip.open(fp, "rt") as fh:
            df = pd.read_csv(fh, sep="\t", low_memory=False)

        # Drop columns we don't need to save RAM
        drop_cols = [c for c in df.columns if c not in _KEEP_COLS]
        df.drop(columns=drop_cols, inplace=True, errors="ignore")

        # Ensure Z exists
        if "Z" not in df.columns and {"BETA","SE"}.issubset(df.columns):
            df["Z"] = df["BETA"].astype(float) / df["SE"].astype(float).clip(lower=1e-12)

        # Downcast numeric columns to float32 — halves memory vs float64
        for col in df.select_dtypes("float64").columns:
            df[col] = df[col].astype("float32")
        for col in df.select_dtypes("int64").columns:
            df[col] = df[col].astype("int32")

        df["cohort"] = cohort
        df["trait"]  = trait
        frames.append(df)
        print(f"    loaded {fp}  →  {len(df):,} rows  trait={trait}")

    combined = pd.concat(frames, ignore_index=True)
    print(f"  Total: {len(combined):,} rows, {combined['trait'].nunique()} traits, "
          f"RAM ≈ {combined.memory_usage(deep=True).sum()/1e9:.2f} GB")
    return combined

def save_tsv(df: pd.DataFrame, path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    print(f"  ✓ saved {path}")

def fig_to_html(fig) -> str:
    return pio.to_html(fig, full_html=False, include_plotlyjs=False,
                       config={"displayModeBar": True, "responsive": True})


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 1 · MULTIVARIATE HYPOTHESIS TESTING                             ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def phase1_hypothesis_testing(ss: pd.DataFrame, out_dir: str) -> dict:
    """
    Phase 1 — Multivariate Hypothesis Testing.

    Memory strategy: we NEVER build a full SNP × trait pivot over all 1.4M SNPs.
    Instead:
      - Correlations  : computed per-cohort by merging pairs of trait series on SNP,
                        then sampling up to MAX_CORR_SNPS shared SNPs.
      - Hotelling T²  : uses only the top-500 genome-wide significant SNPs per cohort.
      - Tornado plots : top-20 loci by |Z| in Cohort A, merged to Cohort B.
    """
    print("\n── Phase 1: Multivariate Hypothesis Testing ──")
    MAX_CORR_SNPS = 200_000   # capforPearson/Spearman
    GWAS_P = 5e-8

    os.makedirs(f"{out_dir}/phase1_hypothesis", exist_ok=True)
    results = {}
    figs    = {}

    cohorts = ss["cohort"].unique().tolist()
    traits  = ss["trait"].unique().tolist()

    # ── 1a. Pairwise_cross-trait_correlations ────
    print("  1a. Cross-trait correlations …")
    corr_results = []
    for cohort in cohorts:
        grp_c = ss[ss["cohort"] == cohort]
        for i, t1 in enumerate(traits):
            for j, t2 in enumerate(traits):
                if j <= i:
                    continue
                s1 = grp_c[grp_c["trait"]==t1][["SNP","Z"]].dropna()
                s2 = grp_c[grp_c["trait"]==t2][["SNP","Z"]].dropna()
                merged = s1.merge(s2, on="SNP", suffixes=("_1","_2"))
                if len(merged) < 30:
                    continue
                # Subsample if very large
                if len(merged) > MAX_CORR_SNPS:
                    merged = merged.sample(MAX_CORR_SNPS, random_state=42)
                x, y = merged["Z_1"].values.astype(float), merged["Z_2"].values.astype(float)
                mask = np.isfinite(x) & np.isfinite(y)
                x, y = x[mask], y[mask]
                n = len(x)
                if n < 30:
                    continue
                r_p, p_p = stats.pearsonr(x, y)
                r_s, p_s = stats.spearmanr(x, y)
                z_r   = np.arctanh(np.clip(r_p, -0.9999, 0.9999))
                se    = 1.0 / np.sqrt(n - 3)
                corr_results.append(dict(
                    cohort=cohort, trait1=t1, trait2=t2,
                    pearson_r=round(float(r_p), 4), pearson_p=float(p_p),
                    spearman_r=round(float(r_s), 4), spearman_p=float(p_s),
                    n_snps=n,
                    ci_lower=round(float(np.tanh(z_r - 1.96*se)), 4),
                    ci_upper=round(float(np.tanh(z_r + 1.96*se)), 4),
                ))
                print(f"    {cohort} {t1} × {t2}  r={r_p:.3f}  n={n:,}")

    corr_df = pd.DataFrame(corr_results)
    save_tsv(corr_df, f"{out_dir}/phase1_hypothesis/cross_trait_correlations.tsv")
    results["cross_trait_corr"] = corr_df

    # ── 1b. Correlation heatmap ───────────────────────────
    print("  1b. Correlation heatmaps …")
    heatmap_figs = {}
    corr_mats    = {}

    for cohort in cohorts:
        sub = corr_df[corr_df["cohort"]==cohort]
        if sub.empty:
            continue
        # Build square matrix from pairwise results
        mat = pd.DataFrame(np.eye(len(traits)), index=traits, columns=traits)
        for _, row in sub.iterrows():
            mat.loc[row.trait1, row.trait2] = row.pearson_r
            mat.loc[row.trait2, row.trait1] = row.pearson_r
        corr_mats[cohort] = mat

        # p-value_matrix
        n_approx = int(sub["n_snps"].median()) if not sub.empty else 50000
        p_mat = pd.DataFrame(np.ones((len(traits), len(traits))),
                              index=traits, columns=traits)
        for _, row in sub.iterrows():
            r = row.pearson_r
            t_s = r * np.sqrt(n_approx-2) / np.sqrt(max(1e-9,1-r**2))
            p_v = 2 * stats.t.sf(abs(t_s), df=n_approx-2)
            p_mat.loc[row.trait1, row.trait2] = p_v
            p_mat.loc[row.trait2, row.trait1] = p_v

        annot = [[
            f"{mat.values[i][j]:.2f}" +
            ("***" if p_mat.values[i][j]<0.001 else
             "**"  if p_mat.values[i][j]<0.01  else
             "*"   if p_mat.values[i][j]<0.05  else "")
            for j in range(len(traits))]
            for i in range(len(traits))]

        fig = go.Figure(go.Heatmap(
            z=mat.values, x=traits, y=traits,
            text=annot, texttemplate="%{text}",
            colorscale=[[0,"#2C5F8A"],[0.5,"#FFFFFF"],[1,"#E05A3A"]],
            zmin=-1, zmax=1,
            colorbar=dict(title="Pearson r", thickness=15),
            hovertemplate="<b>%{y} × %{x}</b><br>r = %{z:.3f}<extra></extra>",
        ))
        apply_base(fig,
            title=dict(text=f"<b>Cross-Trait Genetic Correlation Matrix — {cohort}</b><br>"
                            f"<sup>n ≈ {n_approx:,} shared SNPs · *p<0.05 **p<0.01 ***p<0.001</sup>",
                       x=0.5),
            xaxis=dict(tickangle=-40, showgrid=False),
            yaxis=dict(showgrid=False),
            width=750, height=680)
        heatmap_figs[cohort] = fig_to_html(fig)

    # ── 1c. Hotelling T² — using only genome-wide significant loci ──────────
    print("  1c. Hotelling T² …")
    hotelling_results = []
    for cohort in cohorts:
        grp_c = ss[ss["cohort"]==cohort]
        # Collect top-500 lead SNPs across all traits by |Z|
        sig_snps = (grp_c[grp_c["P"].astype(float) < GWAS_P]
                    .nlargest(500, "Z")[["SNP"]].drop_duplicates())
        if len(sig_snps) < 10:
            # Fall back to top-500 by |Z| if no GW-sig SNPs
            sig_snps = (grp_c.assign(absZ=grp_c["Z"].abs())
                        .nlargest(500, "absZ")[["SNP"]].drop_duplicates())

        # Build small pivot (500 SNPs × n_traits)
        sub_sig = grp_c[grp_c["SNP"].isin(sig_snps["SNP"])]
        pivot = (sub_sig.pivot_table(index="SNP", columns="trait",
                                     values="Z", aggfunc="first")
                        .dropna())
        if pivot.shape[0] < 5 or pivot.shape[1] < 2:
            continue

        # Random background: 500 SNPs sampled from non-significant
        non_sig = grp_c[~grp_c["SNP"].isin(sig_snps["SNP"])]
        rand_snps = non_sig["SNP"].drop_duplicates().sample(
            min(500, non_sig["SNP"].nunique()), random_state=42)
        sub_rand = non_sig[non_sig["SNP"].isin(rand_snps)]
        pivot_r = (sub_rand.pivot_table(index="SNP", columns="trait",
                                        values="Z", aggfunc="first")
                           .dropna())
        shared_traits = list(set(pivot.columns) & set(pivot_r.columns))
        if len(shared_traits) < 2:
            continue

        g1 = pivot[shared_traits].values.astype(float)
        g2 = pivot_r[shared_traits].values.astype(float)
        n1, n2, p = len(g1), len(g2), len(shared_traits)
        S1 = np.cov(g1.T); S2 = np.cov(g2.T)
        Sp = ((n1-1)*S1 + (n2-1)*S2) / (n1+n2-2)
        try:
            Sp_inv = inv(Sp + np.eye(p)*1e-6)
        except Exception:
            continue
        diff = g1.mean(0) - g2.mean(0)
        T2   = n1*n2/(n1+n2) * diff @ Sp_inv @ diff
        F    = T2*(n1+n2-p-1)/(p*(n1+n2-2))
        pval = 1 - stats.f.cdf(F, p, n1+n2-p-1)
        hotelling_results.append(dict(
            cohort=cohort, n_traits=p, n_sig_snps=n1, n_rand_snps=n2,
            T2=round(float(T2),4), F_stat=round(float(F),4),
            df1=p, df2=n1+n2-p-1, p_value=float(pval),
            significant=(pval < 0.05)
        ))
        print(f"    {cohort}: T²={T2:.2f}  F={F:.2f}  p={pval:.3e}")

    hotelling_df = pd.DataFrame(hotelling_results)
    save_tsv(hotelling_df, f"{out_dir}/phase1_hypothesis/hotelling_T2.tsv")
    results["hotelling"] = hotelling_df

    # ── 1d. Cross-cohort tornado plots ──────────────────────────────────────
    print("  1d. Tornado plots …")
    tornado_figs = {}
    qt_traits = [t for t in traits
                 if any(k in t for k in ["ldl","bmi","crp","height"])]
    for t in qt_traits:
        grp_A = ss[(ss["cohort"]=="Cohort_A") & (ss["trait"]==t)]
        grp_B = ss[(ss["cohort"]=="Cohort_B") & (ss["trait"]==t)]
        if grp_A.empty or grp_B.empty:
            continue
        top20 = grp_A.nlargest(20,"Z")[["SNP","CHR","BP","BETA","Z","P"]].copy()
        merged = top20.merge(
            grp_B[["SNP","BETA","Z","P"]].rename(
                columns={"BETA":"BETA_B","Z":"Z_B","P":"P_B"}),
            on="SNP", how="inner")
        merged = merged.rename(columns={"BETA":"BETA_A","Z":"Z_A","P":"P_A"})
        if merged.empty:
            continue
        loci = [f"chr{int(r.CHR)}:{int(r.BP)}" for _,r in merged.iterrows()]
        fig = go.Figure()
        fig.add_trace(go.Bar(name="Cohort A (EUR+EAS)",
            x=merged["BETA_A"].astype(float), y=loci, orientation="h",
            marker_color=PALETTE["cohort_A"],
            hovertemplate="<b>%{y}</b><br>β = %{x:.4f}<extra>Cohort A</extra>"))
        fig.add_trace(go.Bar(name="Cohort B (AFR+AMR+SAS)",
            x=merged["BETA_B"].astype(float), y=loci, orientation="h",
            marker_color=PALETTE["cohort_B"],
            hovertemplate="<b>%{y}</b><br>β = %{x:.4f}<extra>Cohort B</extra>"))
        apply_base(fig, barmode="group",
            title=dict(text=f"<b>Cross-Cohort Effect Size — {t}</b><br>"
                            "<sup>Top 20 loci by |Z| in Cohort A</sup>", x=0.5),
            xaxis=dict(title="Effect size (β)", showgrid=False),
            yaxis=dict(title="Locus", tickfont=dict(size=10), showgrid=False),
            height=600, width=800)
        tornado_figs[t] = fig_to_html(fig)

    figs["heatmaps"]     = heatmap_figs
    figs["tornado"]      = tornado_figs
    figs["hotelling_df"] = hotelling_df
    figs["corr_mats"]    = corr_mats   # passed to static figure generator
    return {"results": results, "figs": figs}


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 2 · REGRESSION ANALYSIS                                         ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def phase2_regression(ss: pd.DataFrame, out_dir: str) -> dict:
    """
    Phase 2 — Regression Analysis.
    Models: OLS, Ridge, Lasso, ElasticNet (α via 5-fold CV), PCR, PLS.
    Target: |Z| (effect-size strength proxy).
    Features: CHR, MAF (if present), −log10(p).
    Each trait×cohort group is capped at MAX_REG_ROWS rows for speed.
    """
    print("\n── Phase 2: Regression Analysis ──")
    MAX_REG_ROWS = 100_000   # more than enough for stable coefficient estimates
    os.makedirs(f"{out_dir}/phase2_regression", exist_ok=True)
    results_all = []
    vif_all     = []
    figs        = {}

    groups = list(ss.groupby(["cohort","trait"]))
    for gi, ((cohort, trait), grp) in enumerate(groups, 1):
        print(f"  [{gi}/{len(groups)}] {cohort} · {trait} …", flush=True)
        grp = grp.dropna(subset=["BETA","SE","P"]).copy()
        if len(grp) < 100:
            continue

        # Subsample to keep runtime manageable
        if len(grp) > MAX_REG_ROWS:
            grp = grp.sample(MAX_REG_ROWS, random_state=42)

        grp["log_p"] = -np.log10(grp["P"].astype(float).clip(1e-300, 1))
        grp["abs_z"] = grp["Z"].abs() if "Z" in grp.columns \
                       else (grp["BETA"].astype(float)/grp["SE"].astype(float).clip(lower=1e-12)).abs()

        feature_cols = ["CHR"]
        if "A1FREQ" in grp.columns:
            grp["maf"] = np.minimum(grp["A1FREQ"].astype(float),
                                    1 - grp["A1FREQ"].astype(float))
            feature_cols.append("maf")
        feature_cols += ["log_p"]

        X_raw = grp[feature_cols].astype(float).values
        y     = grp["abs_z"].astype(float).values
        mask  = np.isfinite(X_raw).all(axis=1) & np.isfinite(y)
        X_raw, y = X_raw[mask], y[mask]
        if len(y) < 50:
            continue

        scaler = StandardScaler()
        X_sc   = scaler.fit_transform(X_raw)
        Xc     = sm.add_constant(X_sc)
        kf     = KFold(n_splits=5, shuffle=True, random_state=42)

        # OLS
        try:
            ols       = sm.OLS(y, Xc).fit()
            ols_r2    = float(ols.rsquared)
            ols_r2adj = float(ols.rsquared_adj)
            ols_fstat = float(ols.fvalue)
            ols_fp    = float(ols.f_pvalue)
        except Exception:
            ols_r2 = ols_r2adj = ols_fstat = ols_fp = np.nan

        # Ridge
        ridge_cv = RidgeCV(alphas=np.logspace(-3, 3, 20), cv=kf)
        ridge_cv.fit(X_sc, y)
        ridge_r2 = float(r2_score(y, ridge_cv.predict(X_sc)))

        # Lasso
        lasso_cv = LassoCV(cv=kf, max_iter=5000, random_state=42, n_alphas=30)
        lasso_cv.fit(X_sc, y)
        lasso_r2 = float(r2_score(y, lasso_cv.predict(X_sc)))

        # ElasticNet
        en_cv = ElasticNetCV(cv=kf, max_iter=5000, random_state=42, n_alphas=20)
        en_cv.fit(X_sc, y)
        en_r2 = float(r2_score(y, en_cv.predict(X_sc)))

        # PCR
        n_comps = min(X_sc.shape[1], 5)
        pca_obj = skPCA(n_components=n_comps)
        X_pca   = pca_obj.fit_transform(X_sc)
        try:
            pcr_fit = sm.OLS(y, sm.add_constant(X_pca)).fit()
            pcr_r2  = float(pcr_fit.rsquared)
        except Exception:
            pcr_r2 = np.nan

        # PLS
        try:
            pls    = PLSRegression(n_components=min(n_comps, 2))
            pls.fit(X_sc, y)
            pls_r2 = float(r2_score(y, pls.predict(X_sc).ravel()))
        except Exception:
            pls_r2 = np.nan

        # VIF
        try:
            for k, fname in enumerate(feature_cols):
                vif_v = variance_inflation_factor(X_sc, k)
                vif_all.append(dict(cohort=cohort, trait=trait,
                                    feature=fname, VIF=round(float(vif_v), 3)))
        except Exception:
            pass

        results_all.append(dict(
            cohort=cohort, trait=trait, n=len(y),
            OLS_R2=round(ols_r2,4), OLS_R2adj=round(ols_r2adj,4),
            OLS_F=round(ols_fstat,3) if not np.isnan(ols_fstat) else np.nan,
            OLS_Fp=ols_fp,
            Ridge_alpha=round(float(ridge_cv.alpha_),4), Ridge_R2=round(ridge_r2,4),
            Lasso_alpha=round(float(lasso_cv.alpha_),4), Lasso_R2=round(lasso_r2,4),
            EN_alpha=round(float(en_cv.alpha_),4),
            EN_l1_ratio=round(float(en_cv.l1_ratio_),4), EN_R2=round(en_r2,4),
            PCR_R2=round(pcr_r2,4), PLS_R2=round(pls_r2,4),
        ))

    reg_df = pd.DataFrame(results_all)
    vif_df = pd.DataFrame(vif_all)
    save_tsv(reg_df, f"{out_dir}/phase2_regression/regression_comparison.tsv")
    save_tsv(vif_df, f"{out_dir}/phase2_regression/vif_diagnostics.tsv")

    # ── Figure: R² comparison grouped bar ───────────────────────────────────
    if not reg_df.empty:
        r2_cols = ["OLS_R2","Ridge_R2","Lasso_R2","EN_R2","PCR_R2","PLS_R2"]
        model_labels = ["OLS","Ridge","Lasso","ElasticNet","PCR","PLS"]
        color_seq = [PALETTE["primary"], PALETTE["secondary"], PALETTE["tertiary"],
                     PALETTE["quaternary"], "#F39C12", "#1ABC9C"]
        fig = go.Figure()
        for col, label, color in zip(r2_cols, model_labels, color_seq):
            if col not in reg_df.columns:
                continue
            fig.add_trace(go.Bar(
                name=label,
                x=[f"{r.cohort}|{r.trait}" for _, r in reg_df.iterrows()],
                y=reg_df[col],
                marker_color=color,
                hovertemplate=f"<b>{label}</b><br>%{{x}}<br>R² = %{{y:.3f}}<extra></extra>",
            ))
        apply_base(fig,
            barmode="group",
            title=dict(text="<b>Model R² Comparison Across Traits & Cohorts</b>",x=0.5),
            xaxis=dict(title="Cohort | Trait", tickangle=-40, showgrid=False),
            yaxis=dict(title="R²", showgrid=False, range=[0, 1]),
            height=520, width=1100,
        )
        figs["r2_comparison"] = fig_to_html(fig)

    # ── Figure: VIF heatmap ──────────────────────────────────────────────────
    if not vif_df.empty:
        vif_pivot = vif_df.groupby(["feature","cohort"])["VIF"].mean().unstack(fill_value=1)
        fig2 = go.Figure(go.Heatmap(
            z=vif_pivot.values,
            x=vif_pivot.columns.tolist(),
            y=vif_pivot.index.tolist(),
            colorscale=[[0,"#FFFFFF"],[0.5,"#F5CBA7"],[1.0,"#E05A3A"]],
            zmin=1, colorbar=dict(title="VIF", thickness=15),
            hovertemplate="<b>%{y}</b> (%{x})<br>VIF = %{z:.2f}<extra></extra>",
        ))
        apply_base(fig2,
            title=dict(text="<b>Variance Inflation Factor (VIF) — Collinearity Diagnostics</b>",x=0.5),
            xaxis=dict(showgrid=False), yaxis=dict(showgrid=False),
            height=400, width=600,
        )
        figs["vif_heatmap"] = fig_to_html(fig2)

    return {"reg_df": reg_df, "vif_df": vif_df, "figs": figs}


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 3 · MULTIPLE-TESTING CORRECTION                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def phase3_multiple_testing(ss: pd.DataFrame, out_dir: str) -> dict:
    """
    For each trait × cohort:
      - Bonferroni, Holm-Bonferroni, BH, BY corrections
      - Storey q-value (via lambda grid)
      - Genomic inflation λ_GC (median chi² / 0.4549)
      - PP-plot (p-p plot) and λ_GC comparison figure
      - Significant-hit counts at each threshold
    """
    print("\n── Phase 3: Multiple-Testing Correction ──")
    os.makedirs(f"{out_dir}/phase3_mtc", exist_ok=True)
    all_mtc  = []
    lambda_gc_records = []
    figs = {}

    for (cohort, trait), grp in ss.groupby(["cohort","trait"]):
        grp = grp.dropna(subset=["P"]).copy()
        p_vals = grp["P"].astype(float).values
        p_vals = p_vals[(p_vals > 0) & (p_vals <= 1)]
        if len(p_vals) < 10:
            continue
        n = len(p_vals)

        # λ_GC
        chi2_obs = stats.chi2.ppf(1 - p_vals, df=1)
        lambda_gc = np.median(chi2_obs) / 0.4549
        lambda_gc_records.append(dict(cohort=cohort, trait=trait,
                                      n_snps=n, lambda_gc=round(lambda_gc, 4)))

        # MTC methods (statsmodels)
        rej_bon, p_bon, _, _  = multipletests(p_vals, method="bonferroni")
        rej_holm, p_holm, _, _ = multipletests(p_vals, method="holm")
        rej_bh, p_bh, _, _    = multipletests(p_vals, method="fdr_bh")
        rej_by, p_by, _, _    = multipletests(p_vals, method="fdr_by")

        # Storey q-value (π0 estimation via lambda grid)
        lambdas = np.arange(0.05, 0.95, 0.05)
        pi0_hat = np.array([np.mean(p_vals >= lam) / (1 - lam) for lam in lambdas])
        pi0_hat = np.clip(pi0_hat, 0, 1)
        pi0     = float(np.median(pi0_hat[-10:])) if len(pi0_hat) >= 10 else 1.0
        pi0     = min(pi0, 1.0)
        ranked  = np.argsort(p_vals)
        q_vals  = np.zeros(n)
        q_vals[ranked[-1]] = pi0 * p_vals[ranked[-1]]
        for idx in reversed(ranked[:-1]):
            q_vals[idx] = min(
                pi0 * n * p_vals[idx] / (np.sum(p_vals <= p_vals[idx])),
                q_vals[ranked[np.where(ranked == idx)[0][0] + 1]]
                if np.where(ranked == idx)[0][0] + 1 < n else 1.0
            )
        n_sig_storey = np.sum(q_vals < 0.05)

        all_mtc.append(dict(
            cohort=cohort, trait=trait, n_snps=n,
            lambda_gc=round(lambda_gc, 4), pi0_storey=round(pi0, 4),
            n_bonferroni=int(rej_bon.sum()),
            n_holm=int(rej_holm.sum()),
            n_BH_FDR=int(rej_bh.sum()),
            n_BY_FDR=int(rej_by.sum()),
            n_storey_q05=int(n_sig_sig if (n_sig_sig := n_sig_storey) else 0),
            threshold_bonferroni=5e-8,
            threshold_BH_min=float(p_bh[rej_bh].max()) if rej_bh.any() else np.nan,
        ))

    mtc_df    = pd.DataFrame(all_mtc)
    lambda_df = pd.DataFrame(lambda_gc_records)
    save_tsv(mtc_df,    f"{out_dir}/phase3_mtc/multiple_testing_correction.tsv")
    save_tsv(lambda_df, f"{out_dir}/phase3_mtc/genomic_inflation.tsv")

    # ── Figure: λ_GC bubble chart ────────────────────────────────────────────
    if not lambda_df.empty:
        fig = go.Figure()
        for cohort in lambda_df["cohort"].unique():
            sub = lambda_df[lambda_df["cohort"]==cohort]
            fig.add_trace(go.Scatter(
                x=sub["trait"], y=sub["lambda_gc"],
                mode="markers+text",
                name=cohort,
                text=[f"{v:.3f}" for v in sub["lambda_gc"]],
                textposition="top center",
                marker=dict(
                    size=sub["n_snps"].apply(lambda x: max(10, min(40, x/500))),
                    color=COHORT_COLORS.get(cohort, PALETTE["neutral"]),
                    opacity=0.85,
                    line=dict(width=1, color="white"),
                ),
                hovertemplate="<b>%{x}</b><br>λ_GC = %{y:.4f}<extra>" + cohort + "</extra>",
            ))
        fig.add_hline(y=1.0, line_dash="dash", line_color="#999999",
                      annotation_text="λ = 1.0 (no inflation)", annotation_position="right")
        apply_base(fig,
            title=dict(text="<b>Genomic Inflation Factor (λ_GC)</b><br>"
                            "<sup>Bubble size ∝ number of tested SNPs</sup>", x=0.5),
            xaxis=dict(title="Trait", showgrid=False),
            yaxis=dict(title="λ_GC", showgrid=False),
            height=480, width=800,
        )
        figs["lambda_gc"] = fig_to_html(fig)

    # ── Figure: significant hits comparison ─────────────────────────────────
    if not mtc_df.empty:
        methods = ["n_bonferroni","n_holm","n_BH_FDR","n_BY_FDR","n_storey_q05"]
        labels  = ["Bonferroni","Holm","BH-FDR","BY-FDR","Storey q<0.05"]
        colors  = [PALETTE["primary"], PALETTE["secondary"], PALETTE["tertiary"],
                   PALETTE["quaternary"],"#F39C12"]
        fig2 = go.Figure()
        for col, label, color in zip(methods, labels, colors):
            if col not in mtc_df.columns:
                continue
            fig2.add_trace(go.Bar(
                name=label,
                x=[f"{r.cohort}|{r.trait}" for _, r in mtc_df.iterrows()],
                y=mtc_df[col],
                marker_color=color,
                hovertemplate=f"<b>{label}</b><br>%{{x}}<br>n hits = %{{y}}<extra></extra>",
            ))
        apply_base(fig2,
            barmode="group",
            title=dict(text="<b>Significant Hits by Multiple-Testing Method</b>",x=0.5),
            xaxis=dict(title="Cohort | Trait", tickangle=-40, showgrid=False),
            yaxis=dict(title="# Significant SNPs", showgrid=False),
            height=500, width=1100,
        )
        figs["sig_hits"] = fig_to_html(fig2)

    # ── Figure: QQ-plot (observed vs expected −log10 p) ─────────────────────
    qq_figs = {}
    for (cohort, trait), grp in ss.groupby(["cohort","trait"]):
        p_vals = grp["P"].dropna().astype(float)
        p_vals = p_vals[(p_vals > 0) & (p_vals <= 1)].values
        if len(p_vals) < 100:
            continue
        p_sort  = np.sort(p_vals)
        n       = len(p_sort)
        exp_p   = (np.arange(1, n+1) - 0.5) / n
        obs_lp  = -np.log10(p_sort)
        exp_lp  = -np.log10(exp_p)
        # Subsample for speed
        stride  = max(1, n // 5000)
        obs_s, exp_s = obs_lp[::stride], exp_lp[::stride]
        lambda_gc_val = mtc_df[(mtc_df.cohort==cohort) & (mtc_df.trait==trait)]["lambda_gc"].values
        lgc_str = f"λ_GC = {lambda_gc_val[0]:.3f}" if len(lambda_gc_val) else ""
        color = COHORT_COLORS.get(cohort, PALETTE["neutral"])
        fig_qq = go.Figure()
        fig_qq.add_trace(go.Scatter(
            x=exp_s, y=obs_s, mode="markers",
            marker=dict(size=3, color=color, opacity=0.5),
            name=f"{cohort} ({lgc_str})",
            hovertemplate="exp: %{x:.2f}<br>obs: %{y:.2f}<extra></extra>",
        ))
        max_val = max(obs_s.max(), exp_s.max()) + 0.5
        fig_qq.add_trace(go.Scatter(
            x=[0, max_val], y=[0, max_val],
            mode="lines", line=dict(color="#999", dash="dash"), name="y=x",
            showlegend=True,
        ))
        apply_base(fig_qq,
            title=dict(text=f"<b>QQ-Plot — {trait}</b><br>"
                            f"<sup>{cohort} · {lgc_str}</sup>", x=0.5),
            xaxis=dict(title="Expected −log₁₀(p)", showgrid=False),
            yaxis=dict(title="Observed −log₁₀(p)", showgrid=False),
            height=480, width=550,
        )
        qq_figs[f"{cohort}_{trait}"] = fig_to_html(fig_qq)

    figs["qq_plots"] = qq_figs
    return {"mtc_df": mtc_df, "lambda_df": lambda_df, "figs": figs}


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 4 · CAUSAL INFERENCE                                            ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def phase4_causal_inference(ss: pd.DataFrame, out_dir: str) -> dict:
    """
    Instrumental Variable (IV) / 2SLS causal inference using GWAS SNPs as
    instruments. Per-trait, per-cohort:
      - Select genome-wide significant SNPs (p < 5e-8) as IVs
      - Compute IV F-statistic (weak-instrument test)
      - Estimate causal effect size via 2SLS proxy
        (ratio estimator when n instruments < 5, IVW otherwise)
      - Heterogeneity: Cochran Q, I²
      - Pleiotropy proxy: MR-Egger intercept
      - Instrument strength: R² per IV
      - Sensitivity: leave-one-out causal estimates
      - Cross-trait causal network (using top instruments)
    """
    print("\n── Phase 4: Causal Inference ──")
    os.makedirs(f"{out_dir}/phase4_causal", exist_ok=True)
    results = []
    loo_results = []
    figs = {}

    GWAS_THRESH = 5e-8

    for (cohort, exposure), grp_exp in ss.groupby(["cohort","trait"]):
        # Use genome-wide significant SNPs as instruments
        ivs = grp_exp[grp_exp["P"] < GWAS_THRESH].copy()
        if len(ivs) < 3:
            continue

        beta_x = ivs["BETA"].astype(float).values
        se_x   = ivs["SE"].astype(float).values

        # For each outcome, attempt IV ratio estimate
        other_traits = [t for t in ss[ss["cohort"]==cohort]["trait"].unique()
                        if t != exposure]
        for outcome in other_traits:
            grp_out = ss[(ss["cohort"]==cohort) & (ss["trait"]==outcome)]
            iv_out  = grp_out[grp_out["SNP"].isin(ivs["SNP"])]
            iv_merged = ivs.merge(
                iv_out[["SNP","BETA","SE"]].rename(
                    columns={"BETA":"BETA_Y","SE":"SE_Y"}),
                on="SNP", how="inner"
            ).dropna(subset=["BETA","SE","BETA_Y","SE_Y"])
            if len(iv_merged) < 3:
                continue

            bx = iv_merged["BETA"].astype(float).values
            se_bx = iv_merged["SE"].astype(float).values
            by = iv_merged["BETA_Y"].astype(float).values
            se_by = iv_merged["SE_Y"].astype(float).values

            # IVW estimate
            weights = 1.0 / (se_by**2 + 1e-12)
            beta_ivw = np.sum(weights * by * bx) / np.sum(weights * bx**2 + 1e-12)
            se_ivw   = np.sqrt(1.0 / (np.sum(weights * bx**2) + 1e-12))
            z_ivw    = beta_ivw / (se_ivw + 1e-12)
            p_ivw    = 2 * stats.norm.sf(abs(z_ivw))

            # Cochran Q (heterogeneity)
            ratio_ests = by / (bx + 1e-12)
            Q = np.sum(weights * (by - beta_ivw * bx)**2)
            df_Q = len(iv_merged) - 1
            p_Q  = 1 - stats.chi2.cdf(Q, df=df_Q)
            I2   = max(0, (Q - df_Q) / (Q + 1e-12)) * 100

            # MR-Egger intercept proxy
            try:
                egger_X = sm.add_constant(bx)
                egger_m = sm.WLS(by, egger_X, weights=weights).fit()
                egger_intercept   = egger_m.params[0]
                egger_intercept_p = egger_m.pvalues[0]
                egger_slope       = egger_m.params[1]
            except Exception:
                egger_intercept = egger_intercept_p = egger_slope = np.nan

            # IV F-stat
            n_approx = 2504
            R2_per_iv = bx**2 / (bx**2 + n_approx * se_bx**2)
            R2_total  = float(R2_per_iv.sum())
            F_iv = (n_approx * R2_total) / (len(iv_merged) * (1 - R2_total + 1e-12))

            results.append(dict(
                cohort=cohort, exposure=exposure, outcome=outcome,
                n_ivs=len(iv_merged),
                beta_IVW=round(beta_ivw, 5), se_IVW=round(se_ivw, 5),
                z_IVW=round(z_ivw, 4), p_IVW=p_ivw,
                CI_lower=round(beta_ivw - 1.96*se_ivw, 5),
                CI_upper=round(beta_ivw + 1.96*se_ivw, 5),
                Q=round(Q, 3), p_Q=round(p_Q, 4), I2=round(I2, 2),
                egger_intercept=round(egger_intercept, 5) if not np.isnan(egger_intercept) else np.nan,
                egger_intercept_p=round(egger_intercept_p, 4) if not np.isnan(egger_intercept_p) else np.nan,
                egger_slope=round(egger_slope, 5) if not np.isnan(egger_slope) else np.nan,
                F_iv=round(F_iv, 2), R2_total=round(R2_total, 5),
                weak_instrument=(F_iv < 10),
            ))

            # Leave-one-out
            for k in range(len(iv_merged)):
                mask = np.ones(len(iv_merged), dtype=bool)
                mask[k] = False
                bx_k  = bx[mask];  by_k   = by[mask]
                se_y_k = se_by[mask]
                w_k   = 1.0 / (se_y_k**2 + 1e-12)
                b_loo = np.sum(w_k * by_k * bx_k) / (np.sum(w_k * bx_k**2) + 1e-12)
                loo_results.append(dict(
                    cohort=cohort, exposure=exposure, outcome=outcome,
                    dropped_snp=iv_merged["SNP"].iloc[k],
                    beta_IVW_LOO=round(b_loo, 5)
                ))

    causal_df = pd.DataFrame(results)
    loo_df    = pd.DataFrame(loo_results)
    save_tsv(causal_df, f"{out_dir}/phase4_causal/iv_causal_estimates.tsv")
    save_tsv(loo_df,    f"{out_dir}/phase4_causal/loo_sensitivity.tsv")

    # ── Figure: Forest plot (IVW estimates with CI) ──────────────────────────
    forest_figs = {}
    if not causal_df.empty:
        for cohort in causal_df["cohort"].unique():
            sub = causal_df[causal_df["cohort"]==cohort].copy()
            sub["label"] = sub["exposure"] + " → " + sub["outcome"]
            sub = sub.dropna(subset=["beta_IVW","CI_lower","CI_upper"])
            sub = sub.sort_values("beta_IVW")
            if sub.empty:
                continue
            colors_f = [PALETTE["secondary"] if p < 0.05 else PALETTE["neutral"]
                        for p in sub["p_IVW"]]
            fig_f = go.Figure()
            for _, row in sub.iterrows():
                color = PALETTE["secondary"] if row.p_IVW < 0.05 else PALETTE["neutral"]
                fig_f.add_trace(go.Scatter(
                    x=[row.CI_lower, row.beta_IVW, row.CI_upper],
                    y=[row.label]*3,
                    mode="lines+markers",
                    line=dict(color=color, width=1.5),
                    marker=dict(size=[4, 10, 4], color=color,
                                symbol=["line-ew","diamond","line-ew"]),
                    name=row.label,
                    showlegend=False,
                    hovertemplate=f"<b>{row.label}</b><br>β={row.beta_IVW:.4f} "
                                  f"[{row.CI_lower:.4f}, {row.CI_upper:.4f}]<br>"
                                  f"p={row.p_IVW:.3e}<extra></extra>",
                ))
            fig_f.add_vline(x=0, line_dash="dash", line_color="#AAAAAA")
            apply_base(fig_f,
                title=dict(text=f"<b>IVW Causal Effect Estimates — {cohort}</b><br>"
                                "<sup>Red = p < 0.05 · Horizontal bars = 95% CI</sup>",
                           x=0.5),
                xaxis=dict(title="Causal effect (β IVW)", showgrid=False),
                yaxis=dict(showgrid=False, tickfont=dict(size=10)),
                height=max(400, 60 * len(sub) + 120),
                width=800,
            )
            forest_figs[cohort] = fig_to_html(fig_f)

    # ── Figure: Causal network heatmap (exposure × outcome) ──────────────────
    network_figs = {}
    if not causal_df.empty:
        for cohort in causal_df["cohort"].unique():
            sub = causal_df[(causal_df["cohort"]==cohort)].copy()
            if sub.empty:
                continue
            all_traits_net = sorted(set(sub["exposure"]) | set(sub["outcome"]))
            net_mat = pd.DataFrame(np.nan, index=all_traits_net, columns=all_traits_net)
            for _, row in sub.iterrows():
                net_mat.loc[row.exposure, row.outcome] = row.beta_IVW
            fig_net = go.Figure(go.Heatmap(
                z=net_mat.values,
                x=net_mat.columns.tolist(),
                y=net_mat.index.tolist(),
                colorscale=[[0,"#2C5F8A"],[0.5,"#FFFFFF"],[1,"#E05A3A"]],
                zmid=0,
                colorbar=dict(title="β IVW", thickness=15),
                hovertemplate="<b>%{y} → %{x}</b><br>β = %{z:.4f}<extra></extra>",
            ))
            apply_base(fig_net,
                title=dict(text=f"<b>Causal Effect Network (IVW) — {cohort}</b>",x=0.5),
                xaxis=dict(title="Outcome", tickangle=-35, showgrid=False),
                yaxis=dict(title="Exposure", showgrid=False),
                height=600, width=700,
            )
            network_figs[cohort] = fig_to_html(fig_net)

    # ── Figure: LOO sensitivity ──────────────────────────────────────────────
    loo_figs = {}
    if not loo_df.empty and not causal_df.empty:
        # Pick one representative pair per cohort
        for cohort in causal_df["cohort"].unique():
            sig_pair = causal_df[
                (causal_df["cohort"]==cohort) &
                (causal_df["p_IVW"] < 0.05)
            ].sort_values("p_IVW")
            if sig_pair.empty:
                continue
            row0 = sig_pair.iloc[0]
            loo_sub = loo_df[
                (loo_df["cohort"]==cohort) &
                (loo_df["exposure"]==row0["exposure"]) &
                (loo_df["outcome"]==row0["outcome"])
            ]
            if loo_sub.empty:
                continue
            fig_loo = go.Figure()
            fig_loo.add_trace(go.Box(
                y=loo_sub["beta_IVW_LOO"],
                boxpoints="all", jitter=0.3,
                marker_color=COHORT_COLORS.get(cohort, PALETTE["neutral"]),
                name="LOO estimates",
                hovertemplate="β_LOO = %{y:.5f}<extra></extra>",
            ))
            fig_loo.add_hline(
                y=row0["beta_IVW"],
                line_dash="dash", line_color=PALETTE["secondary"],
                annotation_text=f"Full IVW β = {row0['beta_IVW']:.4f}",
                annotation_position="right"
            )
            apply_base(fig_loo,
                title=dict(text=f"<b>Leave-One-Out Sensitivity — {row0['exposure']} → {row0['outcome']}</b><br>"
                                f"<sup>{cohort}</sup>", x=0.5),
                yaxis=dict(title="β IVW (LOO)", showgrid=False),
                xaxis=dict(showgrid=False),
                height=460, width=560,
            )
            loo_figs[cohort] = fig_to_html(fig_loo)

    figs["forest"]  = forest_figs
    figs["network"] = network_figs
    figs["loo"]     = loo_figs
    return {"causal_df": causal_df, "loo_df": loo_df, "figs": figs}


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  HTML REPORT ASSEMBLY                                                   ║
# ╚══════════════════════════════════════════════════════════════════════════╝

_PLOTLYJS = None

def get_plotlyjs() -> str:
    """Return the Plotly.js bundle (cached)."""
    global _PLOTLYJS
    if _PLOTLYJS is None:
        _PLOTLYJS = pio.to_html(go.Figure(), full_html=False,
                                include_plotlyjs=True)
        # extract just the script tag
        import re
        m = re.search(r'(<script[^>]*>.*?</script>)', _PLOTLYJS, re.DOTALL)
        _PLOTLYJS = m.group(1) if m else ""
    return _PLOTLYJS


def build_html_report(p1, p2, p3, p4,
                      author: str, institute: str,
                      out_path: str):
    """Assemble all phases into one standalone HTML report."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M")
    plotlyjs = get_plotlyjs()

    def _section(title, icon, content_html):
        return f"""
        <section class="section">
          <h2 class="section-title"><span class="icon">{icon}</span>{title}</h2>
          <div class="section-body">{content_html}</div>
        </section>"""

    def _table(df: pd.DataFrame, max_rows=200) -> str:
        if df is None or df.empty:
            return "<p class='empty'>No results available.</p>"
        df_show = df.head(max_rows).copy()
        # Round floats
        for col in df_show.select_dtypes("float").columns:
            df_show[col] = df_show[col].apply(
                lambda x: f"{x:.3e}" if (abs(x) < 0.001 and x != 0) else f"{x:.4f}"
                if not np.isnan(x) else "—"
            )
        thead = "".join(f"<th>{c}</th>" for c in df_show.columns)
        rows  = ""
        for _, r in df_show.iterrows():
            cells = "".join(f"<td>{v}</td>" for v in r)
            rows += f"<tr>{cells}</tr>"
        ellipsis = (f"<tr><td colspan='{len(df_show.columns)}' style='text-align:center;"
                    f"color:#999'>… {len(df)-max_rows:,} more rows</td></tr>"
                    if len(df) > max_rows else "")
        return (f"<div class='table-wrap'><table class='data-table'>"
                f"<thead><tr>{thead}</tr></thead>"
                f"<tbody>{rows}{ellipsis}</tbody></table></div>")

    def _figs_html(figs_dict: dict) -> str:
        if not figs_dict:
            return ""
        grid = ""
        for key, html in figs_dict.items():
            if html:
                grid += f"<div class='fig-card'>{html}</div>"
        return f"<div class='fig-grid'>{grid}</div>" if grid else ""

    # ── Phase 1 content ──────────────────────────────────────────────────────
    p1f = p1.get("figs", {})
    p1_content = ""
    p1_content += "<h3>Cross-Trait Genetic Correlation (Z-scores)</h3>"
    p1_content += _figs_html(p1f.get("heatmaps", {}))
    p1_content += "<h3>Cross-Cohort Effect Size Comparison</h3>"
    p1_content += _figs_html(p1f.get("tornado", {}))
    p1_content += "<h3>Hotelling T² Test (top loci vs random)</h3>"
    p1_content += _table(p1f.get("hotelling_df", pd.DataFrame()))
    p1_content += "<h3>Pairwise Cross-Trait Correlations</h3>"
    p1_content += _table(p1.get("results", {}).get("cross_trait_corr", pd.DataFrame()))

    # ── Phase 2 content ──────────────────────────────────────────────────────
    p2f = p2.get("figs", {})
    p2_content  = "<h3>Model R² Comparison</h3>"
    p2_content += p2f.get("r2_comparison", "")
    p2_content += "<h3>VIF Collinearity Diagnostics</h3>"
    p2_content += p2f.get("vif_heatmap", "")
    p2_content += "<h3>Regression Results Table</h3>"
    p2_content += _table(p2.get("reg_df", pd.DataFrame()))

    # ── Phase 3 content ──────────────────────────────────────────────────────
    p3f = p3.get("figs", {})
    p3_content  = "<h3>Genomic Inflation (λ_GC)</h3>"
    p3_content += p3f.get("lambda_gc", "")
    p3_content += "<h3>Significant Hits by MTC Method</h3>"
    p3_content += p3f.get("sig_hits", "")
    p3_content += "<h3>QQ-Plots</h3>"
    p3_content += _figs_html(p3f.get("qq_plots", {}))
    p3_content += "<h3>MTC Summary Table</h3>"
    p3_content += _table(p3.get("mtc_df", pd.DataFrame()))

    # ── Phase 4 content ──────────────────────────────────────────────────────
    p4f = p4.get("figs", {})
    p4_content  = "<h3>IVW Causal Effect Forest Plots</h3>"
    p4_content += _figs_html(p4f.get("forest", {}))
    p4_content += "<h3>Causal Effect Network (Heatmap)</h3>"
    p4_content += _figs_html(p4f.get("network", {}))
    p4_content += "<h3>Leave-One-Out Sensitivity</h3>"
    p4_content += _figs_html(p4f.get("loo", {}))
    p4_content += "<h3>Causal Estimates Table</h3>"
    p4_content += _table(p4.get("causal_df", pd.DataFrame()))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
<title>Multivariate GWAS Analysis — CDC 1.0.0</title>
{plotlyjs}
<style>
  /* ── Reset & base ── */
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  html {{ scroll-behavior: smooth; }}
  body {{
    font-family: 'Inter', 'Segoe UI', Arial, sans-serif;
    background: #FFFFFF;
    color: #1A1A2E;
    line-height: 1.6;
    font-size: 14px;
  }}

  /* ── Header ── */
  .page-header {{
    background: linear-gradient(135deg, #1A2940 0%, #2C5F8A 60%, #3A7EC2 100%);
    color: white;
    padding: 48px 60px 40px;
    position: relative;
    overflow: hidden;
  }}
  .page-header::after {{
    content: '';
    position: absolute; bottom: -30px; right: -30px;
    width: 300px; height: 300px;
    border-radius: 50%;
    background: rgba(255,255,255,0.04);
  }}
  .page-header h1 {{
    font-size: 2.2rem; font-weight: 700; letter-spacing: -0.5px;
    margin-bottom: 6px;
  }}
  .page-header .subtitle {{ font-size: 1rem; opacity: 0.8; margin-bottom: 20px; }}
  .meta-row {{
    display: flex; flex-wrap: wrap; gap: 24px;
    font-size: 0.85rem; opacity: 0.75; margin-top: 12px;
  }}
  .meta-pill {{
    background: rgba(255,255,255,0.12);
    border-radius: 20px; padding: 4px 14px;
  }}

  /* ── Navigation tabs ── */
  .nav-tabs {{
    position: sticky; top: 0; z-index: 100;
    background: #FFFFFF;
    border-bottom: 2px solid #E8ECF0;
    display: flex; gap: 0;
    padding: 0 40px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.06);
  }}
  .nav-tab {{
    padding: 14px 22px;
    font-size: 0.9rem; font-weight: 500;
    color: #666; cursor: pointer;
    border-bottom: 3px solid transparent;
    margin-bottom: -2px;
    transition: all 0.2s;
    white-space: nowrap;
  }}
  .nav-tab:hover {{ color: #2C5F8A; }}
  .nav-tab.active {{ color: #2C5F8A; border-bottom-color: #2C5F8A; font-weight: 600; }}

  /* ── Content ── */
  .content {{ padding: 40px 60px 80px; max-width: 1400px; margin: 0 auto; }}

  .tab-pane {{ display: none; }}
  .tab-pane.active {{ display: block; }}

  .section {{ margin-bottom: 48px; }}
  .section-title {{
    font-size: 1.25rem; font-weight: 600; color: #1A1A2E;
    margin-bottom: 20px; padding-bottom: 10px;
    border-bottom: 1px solid #E8ECF0;
    display: flex; align-items: center; gap: 10px;
  }}
  .section-title .icon {{ font-size: 1.4rem; }}
  h3 {{
    font-size: 1rem; font-weight: 600; color: #2C5F8A;
    margin: 24px 0 12px; padding-left: 10px;
    border-left: 3px solid #2C5F8A;
  }}

  /* ── Figure grid ── */
  .fig-grid {{
    display: flex; flex-wrap: wrap; gap: 24px;
    margin: 16px 0;
  }}
  .fig-card {{
    background: #FAFBFC;
    border: 1px solid #E8ECF0;
    border-radius: 10px;
    padding: 12px;
    flex: 1 1 520px;
    min-width: 320px;
  }}

  /* ── Tables ── */
  .table-wrap {{
    overflow-x: auto; margin: 12px 0;
    border: 1px solid #E8ECF0; border-radius: 8px;
  }}
  .data-table {{
    border-collapse: collapse; width: 100%;
    font-size: 12px;
  }}
  .data-table thead {{ background: #F0F4F8; }}
  .data-table th {{
    padding: 10px 14px; text-align: left;
    font-weight: 600; color: #444;
    border-bottom: 2px solid #DDE3EA;
    white-space: nowrap;
  }}
  .data-table td {{
    padding: 8px 14px; border-bottom: 1px solid #F0F4F8;
    font-family: 'JetBrains Mono', monospace;
  }}
  .data-table tr:hover td {{ background: #F7F9FC; }}
  .data-table tr:last-child td {{ border-bottom: none; }}

  /* ── Summary cards ── */
  .summary-cards {{
    display: flex; flex-wrap: wrap; gap: 16px;
    margin: 0 0 36px;
  }}
  .card {{
    flex: 1 1 160px; min-width: 140px;
    background: #F7F9FC;
    border: 1px solid #E8ECF0;
    border-radius: 10px;
    padding: 20px 16px;
    text-align: center;
  }}
  .card .num {{
    font-size: 2rem; font-weight: 700; color: #2C5F8A;
    display: block; line-height: 1.2;
  }}
  .card .lbl {{ font-size: 0.78rem; color: #888; margin-top: 4px; }}

  .empty {{ color: #999; font-style: italic; padding: 20px 0; }}
  .badge {{
    display: inline-block; padding: 2px 8px;
    border-radius: 12px; font-size: 0.75rem;
    background: #E8F4FD; color: #2C5F8A;
    font-weight: 600; margin-left: 8px;
  }}
</style>
</head>
<body>

<!-- ── Header ── -->
<header class="page-header">
  <h1>Multivariate GWAS Statistical Analysis</h1>
  <div class="subtitle">CDC 1.0.0 · Phase A–D Complete Analysis Report</div>
  <div class="meta-row">
    <span class="meta-pill"> {author}</span>
    <span class="meta-pill"> {institute}</span>
    <span class="meta-pill"> Generated: {ts}</span>
    <span class="meta-pill"> Cohort A: EUR+EAS &nbsp;|&nbsp; Cohort B: AFR+AMR+SAS</span>
    <span class="meta-pill"> Traits: LDL, BMI, CRP, Height | CAD, T2D, HTN</span>
  </div>
</header>

<!-- ── Navigation ── -->
<nav class="nav-tabs">
  <div class="nav-tab active" onclick="showTab('overview')"> Overview</div>
  <div class="nav-tab" onclick="showTab('phase1')"> Hypothesis Testing</div>
  <div class="nav-tab" onclick="showTab('phase2')"> Regression</div>
  <div class="nav-tab" onclick="showTab('phase3')"> Multiple Testing</div>
  <div class="nav-tab" onclick="showTab('phase4')"> Causal Inference</div>
</nav>

<!-- ── Content ── -->
<div class="content">

  <!-- Overview tab -->
  <div id="tab-overview" class="tab-pane active">
    <section class="section">
      <h2 class="section-title"><span class="icon"></span>Pipeline Overview</h2>
      <div class="summary-cards">
        <div class="card"><span class="num">2</span><div class="lbl">Cohorts analysed</div></div>
        <div class="card"><span class="num">7</span><div class="lbl">Traits (4 QT + 3 BT)</div></div>
        <div class="card"><span class="num">4</span><div class="lbl">Analysis phases</div></div>
        <div class="card"><span class="num">6</span><div class="lbl">Regression models</div></div>
        <div class="card"><span class="num">5</span><div class="lbl">MTC methods</div></div>
        <div class="card"><span class="num">IV</span><div class="lbl">2SLS causal inference</div></div>
      </div>
      <div class="fig-card" style="max-width:760px">
        <h3 style="border:none;padding:0;color:#1A1A2E;font-size:1rem">Analysis Workflow</h3>
        <table style="width:100%;border-collapse:collapse;font-size:13px;margin-top:12px">
          <thead style="background:#F0F4F8">
            <tr>
              <th style="padding:10px 14px;text-align:left">Phase</th>
              <th style="padding:10px 14px;text-align:left">Methods</th>
              <th style="padding:10px 14px;text-align:left">Outputs</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td style="padding:10px 14px"><strong>Phase 1</strong><br>Hypothesis Testing</td>
              <td style="padding:10px 14px">MANOVA proxy, Hotelling T², Pearson/Spearman cross-trait correlations, Fisher Z CI</td>
              <td style="padding:10px 14px">Correlation heatmaps, cross-cohort forest plots, T² table</td>
            </tr>
            <tr style="background:#FAFBFC">
              <td style="padding:10px 14px"><strong>Phase 2</strong><br>Regression</td>
              <td style="padding:10px 14px">OLS, Ridge, Lasso, ElasticNet (α via 5-fold CV), PCR, PLS, VIF</td>
              <td style="padding:10px 14px">R² comparison chart, VIF heatmap, model table</td>
            </tr>
            <tr>
              <td style="padding:10px 14px"><strong>Phase 3</strong><br>Multiple Testing</td>
              <td style="padding:10px 14px">Bonferroni, Holm, BH-FDR, BY-FDR, Storey q-value, λ_GC</td>
              <td style="padding:10px 14px">QQ-plots, λ_GC bubble chart, significant-hit comparison</td>
            </tr>
            <tr style="background:#FAFBFC">
              <td style="padding:10px 14px"><strong>Phase 4</strong><br>Causal Inference</td>
              <td style="padding:10px 14px">2SLS/IVW, MR-Egger intercept, Cochran Q/I², LOO sensitivity, causal network</td>
              <td style="padding:10px 14px">Forest plots, network heatmap, LOO boxplot, causal table</td>
            </tr>
          </tbody>
        </table>
      </div>
    </section>
  </div>

  <!-- Phase 1 tab -->
  <div id="tab-phase1" class="tab-pane">
    {_section("Phase 1 — Multivariate Hypothesis Testing", "", p1_content)}
  </div>

  <!-- Phase 2 tab -->
  <div id="tab-phase2" class="tab-pane">
    {_section("Phase 2 — Multivariate Regression Analysis", "", p2_content)}
  </div>

  <!-- Phase 3 tab -->
  <div id="tab-phase3" class="tab-pane">
    {_section("Phase 3 — Multiple-Testing Correction", "", p3_content)}
  </div>

  <!-- Phase 4 tab -->
  <div id="tab-phase4" class="tab-pane">
    {_section("Phase 4 — Causal Inference (Mendelian Randomization proxy)", ", p4_content)}
  </div>

</div><!-- /content -->

<script>
function showTab(name) {{
  document.querySelectorAll('.tab-pane').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.nav-tab').forEach(t => t.classList.remove('active'));
  document.getElementById('tab-' + name).classList.add('active');
  event.currentTarget.classList.add('active');
  window.scrollTo({{top: 0, behavior: 'smooth'}});
}}
</script>
</body>
</html>"""

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    print(f"\n✓ Report written → {out_path}")


# ── CLI ─────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="CDC 1.0.0 — Multivariate GWAS Analysis")
    p.add_argument("--sumstats_pattern", required=True,
                   help="Glob pattern for harmonised *_sumstats.tsv.gz files")
    p.add_argument("--out_dir",  required=True, help="Output directory")
    p.add_argument("--author",   default="Nadeem Khan")
    p.add_argument("--institute",default="INRS-Centre Armand-Frappier Santé-Biotechnologie")
    return p.parse_args()


def main():
    args = parse_args()
    print(f"\n{'='*64}")
    print(f"  CDC 1.0.0 — Multivariate GWAS Statistical Analysis")
    print(f"  Author : {args.author}")
    print(f"  Pattern: {args.sumstats_pattern}")
    print(f"  Out    : {args.out_dir}")
    print(f"{'='*64}\n")

    ss = load_sumstats(args.sumstats_pattern)
    print(f"  Traits found: {sorted(ss['trait'].unique())}")

    # ── Phase 1 ─────────────────────────────────────────────────────────────
    p1 = phase1_hypothesis_testing(ss, args.out_dir)
    if _HAS_SF:
        print("  → Phase 1 static figures …")
        corr_mats = p1.get("figs", {}).get("corr_mats", {})
        for cohort, mat in corr_mats.items():
            n_snps = int(p1.get("results",{})
                           .get("cross_trait_corr", pd.DataFrame())
                           .query(f"cohort=='{cohort}'")["n_snps"].median()
                         ) if not p1.get("results",{}).get("cross_trait_corr",
                                                            pd.DataFrame()).empty else 50000
            SF.p1_correlation_heatmap(mat, cohort, n_snps, args.out_dir)
        # Tornado static figures
        ss_a = ss[ss["cohort"]=="Cohort_A"]
        ss_b = ss[ss["cohort"]=="Cohort_B"]
        for t in ss["trait"].unique():
            if not any(k in t for k in ["ldl","bmi","crp","height"]):
                continue
            grp_A = ss_a[ss_a["trait"]==t]
            grp_B = ss_b[ss_b["trait"]==t]
            if grp_A.empty or grp_B.empty:
                continue
            top20 = grp_A.nlargest(20,"Z")[["SNP","CHR","BP","BETA","Z"]].copy()
            merged = top20.merge(
                grp_B[["SNP","BETA"]].rename(columns={"BETA":"BETA_B"}),
                on="SNP", how="inner")
            if merged.empty:
                continue
            merged = merged.rename(columns={"BETA":"BETA_A"})
            merged["locus"] = [f"chr{int(r.CHR)}:{int(r.BP)}"
                               for _, r in merged.iterrows()]
            SF.p1_tornado(merged, t, args.out_dir)
        SF.p1_hotelling_bar(p1["figs"].get("hotelling_df", pd.DataFrame()),
                            args.out_dir)
    if _HAS_ET:
        print("  → Phase 1 Excel tables …")
        ET.write_phase1_excel(
            p1.get("results",{}).get("cross_trait_corr", pd.DataFrame()),
            p1.get("figs",{}).get("hotelling_df", pd.DataFrame()),
            args.out_dir)

    # ── Phase 2 ─────────────────────────────────────────────────────────────
    p2 = phase2_regression(ss, args.out_dir)
    if _HAS_SF:
        print("  → Phase 2 static figures …")
        SF.p2_r2_grouped(p2.get("reg_df", pd.DataFrame()), args.out_dir)
        SF.p2_vif_heatmap(p2.get("vif_df", pd.DataFrame()), args.out_dir)
        SF.p2_coefficient_path(p2.get("reg_df", pd.DataFrame()), args.out_dir)
    if _HAS_ET:
        print("  → Phase 2 Excel tables …")
        ET.write_phase2_excel(
            p2.get("reg_df", pd.DataFrame()),
            p2.get("vif_df", pd.DataFrame()),
            args.out_dir)

    # ── Phase 3 ─────────────────────────────────────────────────────────────
    p3 = phase3_multiple_testing(ss, args.out_dir)
    lambda_df = p3.get("lambda_df", pd.DataFrame())
    if _HAS_SF:
        print("  → Phase 3 static figures …")
        SF.p3_lambda_gc_bar(lambda_df, args.out_dir)
        SF.p3_sig_hits_grouped(p3.get("mtc_df", pd.DataFrame()), args.out_dir)
        SF.p3_pi0_strip(p3.get("mtc_df", pd.DataFrame()), args.out_dir)
        # QQ per trait × cohort
        for (cohort, trait), grp in ss.groupby(["cohort","trait"]):
            p_vals = grp["P"].dropna().astype(float).values
            p_vals = p_vals[(p_vals>0)&(p_vals<=1)]
            if len(p_vals) < 100:
                continue
            lgc_row = lambda_df[
                (lambda_df["cohort"]==cohort)&(lambda_df["trait"]==trait)]
            lgc = float(lgc_row["lambda_gc"].values[0]) if not lgc_row.empty else 1.0
            SF.p3_qq_plot(p_vals, cohort, trait, lgc, args.out_dir)
    if _HAS_ET:
        print("  → Phase 3 Excel tables …")
        ET.write_phase3_excel(
            p3.get("mtc_df", pd.DataFrame()),
            lambda_df, args.out_dir)

    # ── Phase 4 ─────────────────────────────────────────────────────────────
    p4 = phase4_causal_inference(ss, args.out_dir)
    causal_df = p4.get("causal_df", pd.DataFrame())
    loo_df    = p4.get("loo_df",    pd.DataFrame())
    if _HAS_SF:
        print("  → Phase 4 static figures …")
        for cohort in (causal_df["cohort"].unique() if not causal_df.empty else []):
            SF.p4_forest(causal_df, cohort, args.out_dir)
            SF.p4_network_heatmap(causal_df, cohort, args.out_dir)
            SF.p4_loo(loo_df, causal_df, cohort, args.out_dir)
            SF.p4_funnel(causal_df, cohort, args.out_dir)
    if _HAS_ET:
        print("  → Phase 4 Excel tables …")
        ET.write_phase4_excel(causal_df, loo_df, args.out_dir)

    # ── HTML report ──────────────────────────────────────────────────────────
    build_html_report(
        p1, p2, p3, p4,
        author=args.author,
        institute=args.institute,
        out_path=f"{args.out_dir}/Multivariate_Report.html"
    )

    print("\n✓ All phases complete.")
    print(f"  HTML report : {args.out_dir}/Multivariate_Report.html")
    print(f"  Static figs : {args.out_dir}/phase*/static/*.pdf  *.png")
    print(f"  Excel tables: {args.out_dir}/phase*/*.xlsx\n")


if __name__ == "__main__":
    main()

"""
Disease-only sensitivity analysis: 5 PD samples (GDS5646) vs 5 influenza samples (GDS6063).
Same pipeline as make_figures.py but restricted to disease-only subsets.
Outputs: figures/disease_only_results.csv and figures/fig_sensitivity_volcano.png
"""
import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler

os.makedirs("figures", exist_ok=True)
sns.set_theme(style="whitegrid", context="paper")
RED = "#d62728"; BLUE = "#1f77b4"

PD_SAMPLES = ["GSM1318547","GSM1318548","GSM1318549","GSM1318550","GSM1318551"]
FLU_SAMPLES = ["GSM1684096","GSM1684098","GSM1684100","GSM1684102","GSM1684104"]

def load_expr(csv_file, sample_cols):
    df = pd.read_csv(csv_file, usecols=["Gene symbol"] + sample_cols, low_memory=False)
    df["Gene symbol"] = df["Gene symbol"].astype(str).str.upper().str.strip()
    df = df.dropna()
    dedup = df.groupby("Gene symbol")[sample_cols].mean().reset_index()
    return dedup, set(df["Gene symbol"].dropna().unique())

pd_dedup, genes1 = load_expr("GDS5646.csv", PD_SAMPLES)
inf_dedup, genes2 = load_expr("GDS6063.csv", FLU_SAMPLES)
common = genes1 & genes2
print("DISEASE-ONLY (5 PD vs 5 influenza)")
print("shared:", len(common), "PD-unique:", len(genes1-genes2), "flu-unique:", len(genes2-genes1))

# Venn
fig, ax = plt.subplots(figsize=(7.5, 5))
venn2([genes1, genes2], set_labels=["GDS5646 (5 PD)", "GDS6063 (5 influenza)"], ax=ax)
ax.set_title("Disease-only gene-symbol overlap (n=5 each)")
plt.tight_layout(); plt.savefig("figures/fig_sensitivity_venn.png", dpi=200); plt.close()

# Quantile normalize
def quantile_normalize(df):
    sorted_df = pd.DataFrame(np.sort(df.values, axis=0), index=df.index, columns=df.columns)
    mean_per_rank = sorted_df.mean(axis=1)
    mean_per_rank.index = np.arange(1, len(mean_per_rank) + 1)
    ranks = df.rank(method="min").astype(int)
    return ranks.apply(lambda col: col.map(mean_per_rank))

common_idx = pd_dedup.set_index("Gene symbol").index.intersection(inf_dedup.set_index("Gene symbol").index)
pd_common = pd_dedup.set_index("Gene symbol").loc[common_idx, PD_SAMPLES]
inf_common = inf_dedup.set_index("Gene symbol").loc[common_idx, FLU_SAMPLES]
combined = pd.concat([pd_common, inf_common], axis=1)
combined_norm = quantile_normalize(combined)
pd_norm = combined_norm[PD_SAMPLES]; inf_norm = combined_norm[FLU_SAMPLES]

# Welch t-test + BH FDR
common_list = sorted(genes1 & genes2)
results = []
for gene in common_list:
    pv = pd_norm.loc[gene].values if gene in pd_norm.index else None
    iv = inf_norm.loc[gene].values if gene in inf_norm.index else None
    if pv is None or iv is None or len(pv) < 2 or len(iv) < 2:
        continue
    stat, p = ttest_ind(pv, iv, equal_var=False)
    results.append({"Gene": gene, "t_stat": stat, "p_value": p,
                    "mean_PD": float(np.mean(pv)), "mean_Influenza": float(np.mean(iv))})
de_df = pd.DataFrame(results)
de_df["FDR_BH"] = multipletests(de_df["p_value"], method="fdr_bh")[1]
de_df["log2FC"] = np.log2(de_df["mean_PD"] + 1) - np.log2(de_df["mean_Influenza"] + 1)
date_pat = re.compile(r'^\d{1,2}-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)$', re.IGNORECASE)
de_df = de_df[~de_df["Gene"].str.match(date_pat)]
de_df = de_df.sort_values("FDR_BH").reset_index(drop=True)
de_df.to_csv("figures/disease_only_DE_results.csv", index=False)
print("\nTop 12 DE genes (disease-only):\n", de_df.head(12)[["Gene","t_stat","mean_PD","mean_Influenza","FDR_BH","log2FC"]].to_string())

# Correlation (z-score per condition, manuscript method)
sc = StandardScaler()
pd_z = pd.DataFrame(sc.fit_transform(combined_norm[PD_SAMPLES].T), index=PD_SAMPLES, columns=combined_norm[PD_SAMPLES].index).mean()
inf_z = pd.DataFrame(sc.fit_transform(combined_norm[FLU_SAMPLES].T), index=FLU_SAMPLES, columns=combined_norm[FLU_SAMPLES].index).mean()
r = float(np.corrcoef(pd_z, inf_z)[0, 1])
print(f"\nDisease-only Pearson r (z-scored): {r:.4f}")

# Volcano
sig = (de_df["FDR_BH"] < 0.05) & (de_df["FDR_BH"].notna())
pvals = de_df["p_value"].astype(float).clip(lower=1e-300)
fig, ax = plt.subplots(figsize=(7.5, 5.2))
ax.scatter(de_df.loc[~sig, "log2FC"], -np.log10(pvals[~sig]), s=6, c="grey", alpha=0.5, label="ns")
ax.scatter(de_df.loc[sig, "log2FC"], -np.log10(pvals[sig]), s=8, c=RED, alpha=0.7, label="FDR<0.05")
ax.axvline(0, color="k", lw=0.8, alpha=0.4)
try:
    from adjustText import adjust_text
    texts = [ax.text(rr["log2FC"], -np.log10(max(float(rr["p_value"]), 1e-300)), rr["Gene"], fontsize=7) for _, rr in de_df.head(8).iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
except Exception as e:
    print("adjustText skipped:", e)
ax.set_xlabel("log2 fold change (PD vs influenza, disease-only)"); ax.set_ylabel("-log10 p-value")
ax.set_title("Sensitivity analysis: disease-only (5 PD vs 5 influenza)")
ax.legend()
plt.tight_layout(); plt.savefig("figures/fig_sensitivity_volcano.png", dpi=200); plt.close()

print("\nDisease-only sensitivity analysis complete.")

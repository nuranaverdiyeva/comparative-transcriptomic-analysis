"""
Reproduce the manuscript analysis (matching the notebook pipeline) and
generate publication figures. Outputs PNG files to figures/.
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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

os.makedirs("figures", exist_ok=True)
sns.set_theme(style="whitegrid", context="paper")
RED = "#d62728"; BLUE = "#1f77b4"
LAB_PD = "GDS5646\n(PD + healthy)"
LAB_FLU = "GDS6063\n(flu + no-virus)"

# ---- 1. Load & preprocess (notebook cells 5,7,21) ----
def load_expr(csv_file):
    h = pd.read_csv(csv_file, nrows=1)
    gsm = [c for c in h.columns if "GSM" in c]
    df = pd.read_csv(csv_file, usecols=["Gene symbol"] + gsm, low_memory=False)
    df["Gene symbol"] = df["Gene symbol"].str.upper().str.strip()
    df = df.dropna()
    dedup = df.groupby("Gene symbol")[gsm].mean().reset_index()
    return dedup, set(df["Gene symbol"].dropna().unique()), gsm

pd_dedup, genes1, pd_cols = load_expr("GDS5646.csv")
inf_dedup, genes2, inf_cols = load_expr("GDS6063.csv")
common = genes1 & genes2
print("common:", len(common), "unique PD:", len(genes1-genes2), "unique flu:", len(genes2-genes1))

# ---- 2. Venn ----
fig, ax = plt.subplots(figsize=(7.5, 5))
venn2([genes1, genes2], set_labels=["GDS5646 (PD + healthy)", "GDS6063 (influenza + no-virus)"], ax=ax)
ax.set_title("Gene-symbol overlap between datasets")
plt.tight_layout(); plt.savefig("figures/fig1_venn_overlap.png", dpi=200); plt.close()

# ---- 3. Quantile normalize (notebook cell 20) ----
def quantile_normalize(df):
    sorted_df = pd.DataFrame(np.sort(df.values, axis=0), index=df.index, columns=df.columns)
    mean_per_rank = sorted_df.mean(axis=1)
    mean_per_rank.index = np.arange(1, len(mean_per_rank) + 1)
    ranks = df.rank(method="min").astype(int)
    return ranks.apply(lambda col: col.map(mean_per_rank))

common_idx = pd_dedup.set_index("Gene symbol").index.intersection(inf_dedup.set_index("Gene symbol").index)
pd_common = pd_dedup.set_index("Gene symbol").loc[common_idx, pd_cols]
inf_common = inf_dedup.set_index("Gene symbol").loc[common_idx, inf_cols]
combined = pd.concat([pd_common, inf_common], axis=1)
combined_norm = quantile_normalize(combined)
pd_norm = combined_norm[pd_cols]; inf_norm = combined_norm[inf_cols]
print("normalized shape:", combined_norm.shape)

# ---- 4. PCA before/after normalization ----
labels = [LAB_PD] * len(pd_cols) + [LAB_FLU] * len(inf_cols)
def run_pca(df, lab, title, fname):
    vals = df.values.T.astype(float).copy()
    col_mean = np.nanmean(vals, axis=0); col_mean = np.where(np.isnan(col_mean), 0, col_mean)
    idx = np.where(np.isnan(vals)); vals[idx] = np.take(col_mean, idx[1])
    vals = np.nan_to_num(vals, nan=0.0)
    pca = PCA(n_components=2); comps = pca.fit_transform(vals)
    fig, ax = plt.subplots(figsize=(7, 5.2))
    for g, col, lab_txt in [(LAB_PD, RED, "GDS5646"), (LAB_FLU, BLUE, "GDS6063")]:
        ii = [i for i, l in enumerate(lab) if l == g]
        ax.scatter(comps[ii, 0], comps[ii, 1], label=lab_txt, c=col, s=90, edgecolor="k", zorder=3)
    ax.set_xlabel("PC1 (%.1f%%)" % (pca.explained_variance_ratio_[0]*100))
    ax.set_ylabel("PC2 (%.1f%%)" % (pca.explained_variance_ratio_[1]*100))
    ax.set_title(title); ax.legend()
    plt.tight_layout(); plt.savefig(fname, dpi=200); plt.close()

run_pca(combined, labels, "PCA of samples (raw, pre-normalization)", "figures/fig2a_pca_raw.png")
run_pca(combined_norm, labels, "PCA of samples (after joint quantile normalization)", "figures/fig2b_pca_normalized.png")

# ---- 5. Distribution boxplots ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
sns.boxplot(data=combined, ax=axes[0], color="#9ecae1")
axes[0].set_xticklabels([("PD" if l == LAB_PD else "Flu") + "\n" + str(i+1) for i, l in enumerate(labels)], rotation=90, fontsize=6)
axes[0].set_title("Raw expression (shared genes)"); axes[0].set_ylabel("Expression value")
sns.boxplot(data=combined_norm, ax=axes[1], color="#a1d99b")
axes[1].set_xticklabels([("PD" if l == LAB_PD else "Flu") + "\n" + str(i+1) for i, l in enumerate(labels)], rotation=90, fontsize=6)
axes[1].set_title("After joint quantile normalization"); axes[1].set_ylabel("Normalized value")
plt.tight_layout(); plt.savefig("figures/fig3_distribution_boxplots.png", dpi=200); plt.close()

# ---- 6. Welch t-test + BH FDR (notebook cell 24), volcano ----
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
print("Top 15 DE genes:\n", de_df.head(15).to_string())
de_df.to_csv("figures/differential_expression_results.csv", index=False)

sig = (de_df["FDR_BH"] < 0.05) & (de_df["FDR_BH"].notna())
pvals = de_df["p_value"].astype(float).clip(lower=1e-300)
fig, ax = plt.subplots(figsize=(7.5, 5.2))
ax.scatter(de_df.loc[~sig, "log2FC"], -np.log10(pvals[~sig]), s=6, c="grey", alpha=0.5, label="ns")
ax.scatter(de_df.loc[sig, "log2FC"], -np.log10(pvals[sig]), s=8, c=RED, alpha=0.7, label="FDR<0.05")
ax.axvline(0, color="k", lw=0.8, alpha=0.4)
# label only top 8 with adjustText to avoid overlap
try:
    from adjustText import adjust_text
    texts = []
    for _, rr in de_df.head(8).iterrows():
        texts.append(ax.text(rr["log2FC"], -np.log10(max(float(rr["p_value"]), 1e-300)), rr["Gene"], fontsize=7))
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
except Exception as e:
    print("adjustText skipped:", e)
    for _, rr in de_df.head(8).iterrows():
        ax.annotate(rr["Gene"], (rr["log2FC"], -np.log10(max(float(rr["p_value"]), 1e-300))), fontsize=7)
ax.set_xlabel("log2 fold change (GDS5646 vs GDS6063)"); ax.set_ylabel("-log10 p-value")
ax.set_title("Volcano plot: differential expression (Welch t-test + BH FDR)")
ax.legend()
plt.tight_layout(); plt.savefig("figures/fig4_volcano.png", dpi=200); plt.close()

# ---- 7. GO enrichment bar chart (values from manuscript Tables 2 & 3) ----
shared_go = ["Regulation of biological process", "Regulation of cellular process", "Cell communication",
             "Developmental process", "Signaling", "Protein binding"]
flu_go = [("Multicellular organismal process", 109), ("System process", 98),
         ("G protein-coupled receptor signaling", 95), ("Cell periphery", 94), ("Sensory perception", 93)]
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
axes[0].barh(shared_go[::-1], [300]*6, color="#5b9bd5")
axes[0].set_xlabel("-log10 p-value (capped at 300)"); axes[0].set_title("Shared genes - top GO terms")
axes[1].barh([t[0] for t in flu_go][::-1], [t[1] for t in flu_go][::-1], color="#ed7d31")
axes[1].set_xlabel("-log10 p-value"); axes[1].set_title("Influenza-unique genes - top GO terms")
plt.tight_layout(); plt.savefig("figures/fig6_go_enrichment.png", dpi=200); plt.close()

print("\nAll figures saved to figures/")

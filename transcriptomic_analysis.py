"""
Comparative Transcriptomic Analysis of Neurodegeneration and Viral Infection:
A Multi-Method Validation Framework

Compares gene expression profiles from a neurodegenerative disease dataset
(GDS5646, Parkinson's disease) against a viral infection dataset
(GDS6063, Influenza), to identify shared and divergent gene expression
patterns using statistical, correlation, and pathway enrichment methods.
"""

import io
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from bioservices import KEGG
from gprofiler import GProfiler
from matplotlib_venn import venn2
from scipy.stats import chi2_contingency, fisher_exact, shapiro
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

DATA_DIR = "data"
RESULTS_DIR = "results"
BACKGROUND_GENOME_SIZE = 20000  # approximate protein-coding gene count for human


# ============================================================
# Data Extraction & Preprocessing
# ============================================================

def convert_soft_to_csv(soft_file: str, csv_file: str) -> None:
    """Extract the data table from an NCBI GEO SOFT file and save it as CSV."""
    with open(soft_file, "r") as file:
        lines = file.readlines()

    start_idx = None
    for i, line in enumerate(lines):
        if line.startswith("!dataset_table_begin"):
            start_idx = i + 1
            break

    if start_idx is None:
        raise ValueError(f"Table not found in SOFT file: {soft_file}")

    table_data = lines[start_idx:]
    df = pd.read_csv(io.StringIO("".join(table_data)), sep="\t", low_memory=False)
    df.to_csv(csv_file, index=False)
    print(f"Saved: {csv_file}")


def load_expression_data(csv_file: str) -> pd.DataFrame:
    """Load a gene expression CSV, keeping only the gene symbol and GSM sample columns."""
    header = pd.read_csv(csv_file, nrows=1)
    gsm_columns = [col for col in header.columns if "GSM" in col]
    df = pd.read_csv(csv_file, usecols=["Gene symbol"] + gsm_columns)

    df["Gene symbol"] = df["Gene symbol"].str.upper().str.strip()
    return df.dropna()


# ============================================================
# Gene Set Comparison
# ============================================================

def compare_gene_sets(genes_a: set, genes_b: set) -> dict:
    """Compute common and dataset-unique gene sets between two gene symbol sets."""
    common = genes_a & genes_b
    unique_a = genes_a - genes_b
    unique_b = genes_b - genes_a

    print(f"Common genes: {len(common)}")
    print(f"Unique to dataset A: {len(unique_a)}")
    print(f"Unique to dataset B: {len(unique_b)}")

    return {"common": common, "unique_a": unique_a, "unique_b": unique_b}


def save_gene_sets(gene_sets: dict, results_dir: str = RESULTS_DIR) -> None:
    os.makedirs(results_dir, exist_ok=True)
    pd.DataFrame({"Common Genes": sorted(gene_sets["common"])}).to_csv(
        os.path.join(results_dir, "common_genes.csv"), index=False
    )
    pd.DataFrame({"Unique to Dataset A": sorted(gene_sets["unique_a"])}).to_csv(
        os.path.join(results_dir, "unique_dataset_a.csv"), index=False
    )
    pd.DataFrame({"Unique to Dataset B": sorted(gene_sets["unique_b"])}).to_csv(
        os.path.join(results_dir, "unique_dataset_b.csv"), index=False
    )


# ============================================================
# Statistical Analysis
# ============================================================

def fishers_exact_test(common: set, unique_a: set, unique_b: set,
                        background_size: int = BACKGROUND_GENOME_SIZE) -> tuple:
    """
    Two-tailed Fisher's exact test on gene set overlap significance.
    Background size approximates total protein-coding genes not observed in either set.
    """
    a = len(common)
    b = len(unique_a)
    c = len(unique_b)
    d = background_size - (a + b + c)

    table = [[a, b], [c, d]]
    odds_ratio, p_value = fisher_exact(table)

    print(f"Fisher's Exact Test: odds ratio={odds_ratio:.4f}, p={p_value:.6g}")
    print(f"Contingency table: {table}")
    return odds_ratio, p_value


def chi_square_test(table: list) -> float:
    """Chi-square test of independence on a contingency table."""
    chi2, p, _, _ = chi2_contingency(table)
    print(f"Chi-Square Test: p={p:.6g}")
    return p


def shapiro_normality_test(values: pd.Series) -> float:
    """Shapiro-Wilk test for normality of a single expression column."""
    stat, p = shapiro(values)
    print(f"Shapiro-Wilk Test: statistic={stat:.4f}, p={p:.6g}")
    return p


def pearson_correlation_matrix(expression_data: pd.DataFrame) -> pd.DataFrame:
    """Gene-wise Pearson correlation matrix across samples."""
    return expression_data.T.corr(method="pearson")


# ============================================================
# Cross-Condition Comparison
# ============================================================

def normalize_and_compare_conditions(df_a: pd.DataFrame, df_b: pd.DataFrame,
                                      common_genes: list) -> pd.DataFrame:
    """
    Z-score normalize both conditions' expression data (per sample), align on
    shared genes, and return a merged DataFrame of per-gene mean expression
    for each condition.
    """
    df_a = df_a.set_index("gene symbol") if "gene symbol" in df_a.columns else df_a
    df_b = df_b.set_index("gene symbol") if "gene symbol" in df_b.columns else df_b

    df_a = df_a.groupby(df_a.index).mean(numeric_only=True)
    df_b = df_b.groupby(df_b.index).mean(numeric_only=True)

    shared = sorted(set(df_a.index).intersection(df_b.index).intersection(common_genes))
    expr_a = df_a.loc[shared]
    expr_b = df_b.loc[shared]

    scaler = StandardScaler()
    scaled_a = pd.DataFrame(scaler.fit_transform(expr_a.T), index=expr_a.columns, columns=expr_a.index)
    scaled_b = pd.DataFrame(scaler.fit_transform(expr_b.T), index=expr_b.columns, columns=expr_b.index)

    merged = pd.concat([
        scaled_a.mean().rename("Condition_A"),
        scaled_b.mean().rename("Condition_B"),
    ], axis=1)

    correlation = merged.corr().iloc[0, 1]
    print(f"Pearson correlation between conditions (across shared genes): {correlation:.3f}")

    merged["abs_diff"] = (merged["Condition_A"] - merged["Condition_B"]).abs()
    return merged


# ============================================================
# GO Enrichment & KEGG Pathway Analysis
# ============================================================

def run_go_enrichment(gene_list: list, organism: str = "hsapiens") -> pd.DataFrame:
    """Run Gene Ontology enrichment analysis via g:Profiler."""
    gp = GProfiler(return_dataframe=True)
    return gp.profile(organism=organism, query=list(gene_list))


def run_kegg_pathway_analysis(gene_list: list, limit: int = 100) -> pd.DataFrame:
    """
    Look up KEGG pathways for each gene in gene_list (capped at `limit` genes,
    since each gene requires two sequential KEGG API calls).
    """
    kegg = KEGG()
    found_pathways = []

    for gene_symbol in tqdm(gene_list[:limit]):
        try:
            gene_info = kegg.find("genes", gene_symbol)
            if not gene_info:
                continue
            for line in gene_info.strip().split("\n"):
                if line.startswith("hsa:"):
                    kegg_id = line.split("\t")[0]
                    parsed = kegg.parse(kegg.get(kegg_id))
                    if "PATHWAY" in parsed:
                        for path_id, path_name in parsed["PATHWAY"].items():
                            found_pathways.append({
                                "Gene": gene_symbol,
                                "KEGG_ID": kegg_id,
                                "Pathway_ID": path_id,
                                "Pathway_Name": path_name,
                            })
                    break
        except Exception as e:
            print(f"Error with gene {gene_symbol}: {e}")

    return pd.DataFrame(found_pathways)


def get_pathways_for_gene(gene_symbol: str) -> dict:
    """
    Look up KEGG pathways for a single gene of interest.
    Standalone utility, independent of the main gene-set pipeline.
    """
    kegg = KEGG()
    gene_info = kegg.find("genes", gene_symbol)
    if not gene_info:
        print(f"Gene {gene_symbol} not found in KEGG.")
        return {}

    kegg_id = None
    for line in gene_info.split("\n"):
        if line.startswith("hsa:"):
            kegg_id = line.split("\t")[0]
            break

    if not kegg_id:
        print(f"Could not identify a KEGG ID for {gene_symbol}.")
        return {}

    parsed = kegg.parse(kegg.get(kegg_id))
    pathways = parsed.get("PATHWAY", {})

    if pathways:
        print(f"Pathways for {gene_symbol}:")
        for path_id, path_name in pathways.items():
            print(f"  - {path_id}: {path_name}")
    else:
        print(f"No pathways found for {gene_symbol}.")

    return pathways


# ============================================================
# Visualization
# ============================================================

def plot_venn_diagram(genes_a: set, genes_b: set, labels: tuple = ("Dataset A", "Dataset B")) -> None:
    plt.figure(figsize=(10, 6))
    venn2([genes_a, genes_b], set_labels=labels)
    plt.title("Gene Overlap Between Datasets")
    plt.show()


def plot_expression_heatmap(merged_expr: pd.DataFrame, title: str) -> None:
    plt.figure(figsize=(12, 4))
    sns.heatmap(merged_expr.drop(columns="abs_diff", errors="ignore").T,
                cmap="vlag", annot=False, cbar=True)
    plt.title(title)
    plt.xlabel("Genes")
    plt.ylabel("Condition")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def plot_top_differential_genes(merged_expr: pd.DataFrame, top_n: int = 20) -> None:
    top_diff = merged_expr.sort_values("abs_diff", ascending=False).head(top_n).drop(columns="abs_diff")
    plt.figure(figsize=(12, 4))
    sns.heatmap(top_diff.T, cmap="coolwarm", annot=False, cbar=True)
    plt.title(f"Top {top_n} Differentiated Common Genes")
    plt.xlabel("Genes")
    plt.ylabel("Condition")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()


def plot_go_enrichment(go_results: pd.DataFrame, top_n: int = 10) -> None:
    top_pathways = go_results.nlargest(top_n, "p_value")
    plt.figure(figsize=(8, 5))
    sns.barplot(y=top_pathways["name"], x=-np.log10(top_pathways["p_value"]), color="steelblue")
    plt.xlabel("-log10(p-value)")
    plt.ylabel("GO Term")
    plt.title(f"Top {top_n} Enriched GO Terms")
    plt.tight_layout()
    plt.show()


def plot_kegg_pathway_counts(kegg_results: pd.DataFrame) -> None:
    counts = kegg_results.groupby("Gene").size()
    count_per_pathway_bucket = counts.value_counts().sort_index()

    plt.figure(figsize=(8, 5))
    count_per_pathway_bucket.plot(kind="bar", color="mediumorchid", edgecolor="black")
    plt.title("Number of KEGG Pathways per Gene")
    plt.xlabel("KEGG Pathways per Gene")
    plt.ylabel("Number of Genes")
    plt.tight_layout()
    plt.show()


# ============================================================
# Main Pipeline
# ============================================================

def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # --- 1. Extract & preprocess ---
    convert_soft_to_csv(os.path.join(DATA_DIR, "GDS5646_full.soft"), os.path.join(DATA_DIR, "GDS5646.csv"))
    convert_soft_to_csv(os.path.join(DATA_DIR, "GDS6063_full.soft"), os.path.join(DATA_DIR, "GDS6063.csv"))

    parkinson_df = load_expression_data(os.path.join(DATA_DIR, "GDS5646.csv"))
    influenza_df = load_expression_data(os.path.join(DATA_DIR, "GDS6063.csv"))

    # --- 2. Compare gene sets ---
    genes_a = set(parkinson_df["Gene symbol"].unique())
    genes_b = set(influenza_df["Gene symbol"].unique())
    gene_sets = compare_gene_sets(genes_a, genes_b)
    save_gene_sets(gene_sets)
    plot_venn_diagram(genes_a, genes_b, labels=("Parkinson's (GDS5646)", "Influenza (GDS6063)"))

    common_genes = sorted(gene_sets["common"])

    # --- 3. Statistical validation of overlap significance ---
    fishers_exact_test(gene_sets["common"], gene_sets["unique_a"], gene_sets["unique_b"])

    contingency_table = [
        [len(gene_sets["common"]), len(gene_sets["unique_a"])],
        [len(gene_sets["unique_b"]), BACKGROUND_GENOME_SIZE],
    ]
    chi_square_test(contingency_table)

    gsm_columns = [c for c in parkinson_df.columns if c.startswith("GSM")]
    if gsm_columns:
        shapiro_normality_test(parkinson_df[gsm_columns[0]])

    # --- 4. Correlation & cross-condition comparison ---
    merged_expr = normalize_and_compare_conditions(parkinson_df.rename(columns=str.lower),
                                                     influenza_df.rename(columns=str.lower),
                                                     common_genes)
    plot_expression_heatmap(merged_expr, "Common Gene Expression (Mean): Parkinson's vs Influenza")
    plot_top_differential_genes(merged_expr)

    # --- 5. GO enrichment ---
    go_results = run_go_enrichment(common_genes)
    go_results.to_csv(os.path.join(RESULTS_DIR, "GO_enrichment_common_genes.csv"), index=False)
    plot_go_enrichment(go_results)

    # --- 6. KEGG pathway analysis ---
    kegg_results = run_kegg_pathway_analysis(common_genes, limit=100)
    kegg_results.to_csv(os.path.join(RESULTS_DIR, "KEGG_pathway_analysis_common_genes.csv"), index=False)
    if not kegg_results.empty:
        plot_kegg_pathway_counts(kegg_results)

    print("\nPipeline complete. Results saved to the 'results/' directory.")


if __name__ == "__main__":
    main()

# Comparative Transcriptomic Analysis of Neurodegeneration and Viral Infection: A Multi-Method Validation Framework

![Python](https://img.shields.io/badge/python-3.9-blue) ![License](https://img.shields.io/badge/license-MIT-green) ![Status](https://img.shields.io/badge/status-active-brightgreen) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21851016.svg)](https://doi.org/10.5281/zenodo.21851016)

A gene expression analysis pipeline comparing a neurodegenerative disease dataset against a viral infection dataset, using cross-dataset normalization, differential expression testing, correlation analysis, statistical validation, and pathway enrichment to identify shared and divergent transcriptomic signatures.

International collaboration with **Leyla Baghirzada, MD, FRCPC, MPH**, Clinical Assistant Professor, Department of Anesthesiology, Perioperative and Pain Medicine, University of Calgary.

Manuscript submitted to Cifra. Computer Sciences and Informatics (under review). Preprint available on Zenodo: https://doi.org/10.5281/zenodo.21851016

> **Revision (reviewer corrections).** The repository has been updated to address the reviewer comments: a dataset-characteristics table, a batch–condition-confound discussion, deeper GO interpretation, new figures, and direct script permalinks. See [Revision notes](#revision-notes) below.

## Direct links to scripts and figures

- Analysis notebook: [`transcriptomic_analysis.ipynb`](transcriptomic_analysis.ipynb)
- Standalone script: [`transcriptomic_analysis.py`](transcriptomic_analysis.py)
- Figure-generation script: [`make_figures.py`](make_figures.py)
- Sensitivity-analysis script: [`make_figures_sensitivity.py`](make_figures_sensitivity.py)
- Generated figures: [`figures/`](figures)
- Full differential-expression results: [`figures/differential_expression_results.csv`](figures/differential_expression_results.csv)
- Disease-only sensitivity results: [`figures/disease_only_DE_results.csv`](figures/disease_only_DE_results.csv)

## Table of Contents
- [Overview](#overview)
- [Project Structure](#project-structure)
- [Data](#data)
- [Methodology](#methodology)
- [Setup](#setup)
- [Usage](#usage)
- [Key Findings](#key-findings)
- [Limitations](#limitations)
- [Tech Stack](#tech-stack)
- [License](#license)
- [Authors](#authors)

## Overview

This pipeline compares gene expression profiles from Parkinson's disease (GDS5646) and Influenza (GDS6063) datasets, sourced from NCBI's Gene Expression Omnibus, to characterize the overlap between neurodegenerative and viral-response gene expression signatures. The framework applies multiple independent validation methods, statistical significance testing, formal differential expression analysis, correlation analysis, and functional pathway enrichment, so that any observed overlap is supported by more than one line of evidence rather than a single test.

## Project Structure

```
comparative-transcriptomic-analysis/
├── transcriptomic_analysis.ipynb # Current, corrected end-to-end analysis notebook
├── transcriptomic_analysis.py    # Equivalent standalone Python script
├── make_figures.py              # Reproduces all manuscript figures into figures/
├── make_figures_sensitivity.py # Disease-only (5 PD vs 5 influenza) sensitivity analysis
├── figures/                     # Generated figures (PNG) + differential_expression_results.csv + disease_only_DE_results.csv
├── GDS5646_full.soft             # Raw GEO SOFT file — Parkinson's disease dataset
├── GDS6063_full.soft             # Raw GEO SOFT file — Influenza dataset
├── GDS5646.csv                   # Parsed, tabular version of GDS5646
├── GDS6063.csv                   # Parsed, tabular version of GDS6063
├── requirements.txt
├── LICENSE.md
└── README.md
```

> Note: a `.github/workflows/` CI workflow is referenced in earlier versions; if absent, it can be ignored.

## Data

Two public datasets from NCBI's Gene Expression Omnibus. **Both DataSets were generated on the same microarray platform (GPL10558, Illumina HumanHT-12 V4.0 expression beadchip)**; the principal technical difference is sample source/tissue and study-level handling, not array hardware. Each DataSet also contains both case and control samples (see *Group composition* below).

### Table 1. Dataset characteristics

| Feature | GDS5646 (Parkinson's) | GDS6063 (Influenza) |
|---|---|---|
| GEO DataSet accession | [GDS5646](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5646) | [GDS6063](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS6063) |
| Reference series | GSE54536 | GSE68849 |
| Platform (GPL) | GPL10558 | GPL10558 |
| Platform name | Illumina HumanHT-12 V4.0 beadchip | Illumina HumanHT-12 V4.0 beadchip |
| Sample source / tissue | Peripheral blood (whole blood) | Plasmacytoid dendritic cells (pDC), ex vivo |
| Organism | Homo sapiens | Homo sapiens |
| Samples (total) | 10 | 10 |
| Group composition | 5 stage-1 PD + 5 healthy controls | 5 influenza A-exposed + 5 no-virus controls |
| PubMed ID | 24804238 | 26826244 |
| Publication year | 2014 | 2016 |

| Dataset | Condition | Samples | Role |
|---|---|---|---|
| [GDS5646](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5646) | Parkinson's disease | 10 | Neurodegeneration proxy |
| [GDS6063](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS6063) | Influenza infection | 10 | Viral infection proxy |

**Data files included in this repository:**
- `GDS5646_full.soft`, `GDS6063_full.soft` — Raw GEO SOFT-format files downloaded directly from NCBI Gene Expression Omnibus for each dataset.
- `GDS5646.csv`, `GDS6063.csv` — Parsed, tabular versions of the above, generated by `convert_soft_to_csv()` in `transcriptomic_analysis.ipynb`, restricted to GSM-prefixed sample columns and the Gene symbol column, used as the input for all downstream analysis.

## Methodology

**1. Preprocessing**
Parses each SOFT file's data table, standardizes gene symbols (uppercase, whitespace-stripped), drops incomplete records, and collapses duplicate gene symbols (multiple probes mapping to the same gene) by averaging. Gene symbols corrupted into calendar-date strings by prior spreadsheet handling (a known artifact affecting genes such as MARCH1 and the SEPT family) are identified and excluded.

**2. Gene Set Comparison**
Identifies genes common to both conditions versus genes unique to each, as the starting point for downstream analysis.

**3. Cross-Dataset Normalization**
Although GDS5646 and GDS6063 use the same platform (GPL10558), they originate from independent studies, tissues (peripheral blood vs plasmacytoid dendritic cells), and processing environments. Expression values for genes common to both datasets are combined into a single matrix across all 20 samples and jointly quantile-normalized, so both datasets are brought onto a shared expression-value distribution before any cross-condition comparison. Formal batch correction (e.g., ComBat) was **not** applied as a primary inferential step, because dataset of origin and biological condition are perfectly confounded (all PD samples from GDS5646, all influenza samples from GDS6063); see the manuscript's *Technical comparability and the batch–condition confound* subsection.

**4. Statistical Validation**
- **Fisher's Exact Test:** tests whether the observed gene overlap is significant against a constructed reference background, and separately against the observed union of both gene sets
- **Chi-Square Test:** independence test on the same contingency structure, as a cross-check against the Fisher's exact result
- **Shapiro-Wilk Test:** checks normality of expression values, informing the choice of non-parametric and rank-based downstream methods

**5. Differential Expression Analysis**
A per-gene Welch's t-test (unequal variance) compares the 10 GDS5646 samples (5 Parkinson's + 5 healthy controls) against the 10 GDS6063 samples (5 influenza-exposed + 5 no-virus controls) across the jointly-normalized, shared gene set. Resulting p-values are corrected for multiple testing using the Benjamini-Hochberg false discovery rate (FDR) procedure, with log2 fold change computed per gene.

**6. Correlation Analysis**
Expression values are z-score normalized per condition (via `StandardScaler`) before computing Pearson correlation between conditions across shared genes, to assess whether gene-level overlap corresponds to coordinated expression.

**7. Functional Enrichment (GO)**
Runs Gene Ontology enrichment separately on the common gene set and on the influenza-unique gene set via g:Profiler, to characterize the biological processes each set is involved in.

**8. Pathway Analysis (KEGG)**
Cross-references a sample of common genes against KEGG pathway data via the KEGG REST API (through `bioservices`) to identify shared biological pathways. Capped per run (default 100 genes) since each lookup requires sequential API calls. *Note: the KEGG REST API has experienced intermittent outages (HTTP 400 errors); some pathway results in the current analysis are drawn from an earlier successful query window rather than the most recent run.*

## Setup

```bash
pip install -r requirements.txt
```

Place the raw GEO SOFT files in the repository root (already included in this repo):
```
GDS5646_full.soft
GDS6063_full.soft
```

## Usage

Open and run `transcriptomic_analysis.ipynb` top to bottom (Restart Kernel + Run All is recommended to avoid stale variable state between cells). This regenerates:
- `GDS5646.csv`, `GDS6063.csv` — parsed expression tables
- `common_genes.csv`, `unique_GDS5646.csv`, `unique_GDS6063.csv` — gene set comparisons
- `differential_expression_results.csv` — full differential expression results
- `GO_enrichment_common_genes.csv`, `GO_enrichment_GDS6063.csv` — functional enrichment results
- `KEGG_pathway_analysis_common_genes_TEST.csv` — pathway mapping (when the KEGG API is available)

To regenerate the manuscript figures (Venn, PCA before/after normalization, distribution boxplots, volcano plot, GO enrichment):
```bash
pip install -r requirements.txt
python make_figures.py
```
This writes the PNGs and `differential_expression_results.csv` into `figures/`.

To run the disease-only sensitivity analysis (5 PD vs 5 influenza, Supplementary Table S1 / Figures S1–S2):
```bash
python make_figures_sensitivity.py
```
This writes `figures/fig_sensitivity_venn.png`, `figures/fig_sensitivity_volcano.png`, and `figures/disease_only_DE_results.csv`. The sensitivity analysis confirms the two core findings (gene-overlap containment; no coordinated expression, r ≈ 0.005) and that all highlighted immune-signaling genes (RNF125, CD52, HCST, ITPKB, CRBN) remain significantly differentially expressed (FDR < 0.001).

## Key Findings

- **12,579 genes** shared between the two datasets, significantly more than expected by chance (Fisher's exact and chi-square tests, both p < 0.001)
- **Low correlation (r = 0.009)** between conditions across shared genes despite the large overlap, indicating no coordinated expression relationship
- **Differential expression analysis** identified genes significantly distinguishing the two conditions, including several with established roles in antiviral and immune signaling (RNF125, CD52, HCST, ITPKB)
- **Functional enrichment** of the influenza-unique gene set showed system process, sensory perception, and G protein-coupled receptor signaling terms — findings discussed in the context of known intersections between antiviral innate immune pathways and Parkinson's disease risk genes in the accompanying manuscript

## Limitations

- Each dataset serves as a proxy for a broader condition category (a single Parkinson's dataset for "neurodegeneration," a single Influenza dataset for "viral infection"); findings reflect these specific datasets rather than the categories at large. Each DataSet also pools case and control samples (5 PD + 5 healthy; 5 influenza + 5 no-virus), which is disclosed in the dataset-characteristics table and limitations
- No matched healthy-control samples are included within the cross-condition comparison; differential expression compares the two pooled conditions directly rather than each against controls
- Formal batch-correction methods (e.g., ComBat) were not applied as a primary inferential step because dataset of origin and biological condition are perfectly confounded (all PD samples from GDS5646, all influenza from GDS6063); joint quantile normalization was used instead to align marginal distributions
- Sample size (n=10 per condition; effectively n=5 per disease subgroup) limits statistical power for individual-gene-level claims
- KEGG pathway analysis is capped at a configurable gene limit per run and subject to third-party API availability
- No automated test suite currently covers the preprocessing or statistical functions

## Tech Stack

| Category | Tools |
|---|---|
| Data processing | pandas, numpy |
| Statistics | scipy (Fisher's exact, chi-square, Shapiro-Wilk, Welch's t-test), statsmodels (FDR correction) |
| Machine learning utilities | scikit-learn (StandardScaler) |
| Pathway & enrichment analysis | g:Profiler (gprofiler-official), bioservices (KEGG) |
| Visualization | matplotlib, seaborn, matplotlib-venn |

## License

MIT, see [LICENSE.md](LICENSE.md).

## Authors

Nurana Verdiyeva, Istanbul Technical University — in collaboration with Leyla Baghirzada, MD, FRCPC, MPH, Clinical Assistant Professor, Department of Anesthesiology, Perioperative and Pain Medicine, University of Calgary.

## Revision notes (reviewer corrections)

This revision addresses every point raised in the Cifra reviewer report:

1. **Sample size.** Limitations now state n = 10 per DataSet (effectively n = 5 per disease subgroup) and that top DE genes are candidate findings, not biomarkers.
2. **Batch correction / ComBat.** A new Methods subsection (*Technical comparability and the batch–condition confound*) explains that ComBat was deliberately not applied as a primary inferential step because dataset of origin and biological condition are perfectly confounded; joint quantile normalization aligns marginal distributions, and a PCA QC (Figure 2) is provided. Differential-expression results are framed as suggestive, not batch-corrected.
3. **GO interpretation.** The influenza-unique *sensory perception* / *G protein-coupled receptor signaling* enrichment is now discussed as potentially platform-coverage- and tissue-driven (plasmacytoid dendritic cells), and is not advanced as evidence of a shared neuroimmune mechanism.
4. **Biological vs technical split.** The Discussion now explicitly separates candidate biological signal (RNF125, CD52, HCST, ITPKB) from technical/annotation artifacts (near-complete PD-set containment in influenza set; broad shared-gene GO terms; influenza-unique receptor-family enrichment).
5. **Dataset-characteristics table.** A new Table 1 reports platform, tissue, group composition, PubMed ID, and publication year for each DataSet — and documents that both DataSets use the **same** platform (GPL10558).
6. **Direct script links.** The manuscript's *Data and Code Availability* statement and this README now give full permalinks to `transcriptomic_analysis.ipynb`, `transcriptomic_analysis.py`, `make_figures.py`, and `figures/`, plus the Zenodo DOI.

Additional corrections made during revision: the prior "different microarray platforms" / "cross-platform normalization" wording has been corrected to "same platform (GPL10558), different studies and tissues" / "cross-dataset normalization"; table numbering has been reordered to first-appearance order; the "fifteen genes vs eleven rows" inconsistency in the DE table has been flagged for resolution; and the KEGG reproducibility caveat has been softened. See `reviewer_response_and_changes_summary.md` (in the manuscript revision package) for the full point-by-point response.

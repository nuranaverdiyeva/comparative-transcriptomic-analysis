# Revised Manuscript — Changes to Paste into the Cifra Editor

**Article:** Comparative Transcriptomic Analysis of Neurodegeneration and Viral Infection: A Multi-Method Validation Framework
**Authors:** Nurana Verdiyeva¹; Leyla Baghirzada, MD, FRCPC, MPH²

This document contains the **exact replacement/insertion text** for each section that must change to address the reviewer's comments. Each block is labelled `REPLACE …` or `INSERT …` so you can paste it directly into the Cifra editor. New figures are listed at the end with the filenames to upload.

---

## A. NEW — Table 1. Characteristics of the source datasets

> INSERT as a new table at the end of the **Data sources** subsection in Methods (immediately after the sentence ending "...converted to tabular CSV format by parsing the embedded data table.").

**Table 1. Characteristics of the two NCBI GEO DataSets used in this study.**

| Feature | GDS5646 (Parkinson's disease) | GDS6063 (Influenza infection) |
|---|---|---|
| GEO DataSet accession | GDS5646 | GDS6063 |
| Reference series | GSE54536 | GSE68849 |
| Platform (GPL accession) | GPL10558 | GPL10558 |
| Platform name | Illumina HumanHT-12 V4.0 expression beadchip | Illumina HumanHT-12 V4.0 expression beadchip |
| Technology | oligonucleotide beads | oligonucleotide beads |
| Sample source / tissue | Peripheral blood (whole blood) | Primary plasmacytoid dendritic cells (pDC), ex vivo |
| Organism | Homo sapiens | Homo sapiens |
| Sample type | RNA | RNA |
| Total samples | 10 | 10 |
| Group composition within the DataSet | 5 stage-1 PD patients + 5 paired healthy controls | 5 influenza A-exposed + 5 no-virus controls |
| PubMed ID | 24804238 | 26826244 |
| Publication year | 2014 | 2016 |

*Table note.* Although the two DataSets were curated from independent studies, tissues, and processing environments, they were **both generated on the same microarray platform (GPL10558)**. The principal technical difference is therefore sample/tissue origin and study-level handling rather than array hardware. This table makes the cross-study, cross-tissue nature of the comparison explicit before interpreting downstream results.

---

## B. REPLACE the "Cross-platform normalization" subsection (Methods)

> REPLACE the existing "Cross-platform normalization" paragraph with the following two paragraphs.

**Cross-dataset normalization.** Although GDS5646 and GDS6063 were generated on the same microarray platform (GPL10558; Table 1), they originate from independent studies, tissues (peripheral blood versus plasmacytoid dendritic cells), and processing environments, so raw expression values are not directly comparable across the two datasets. To reduce these distributional differences, expression values for the genes common to both datasets were combined into a single matrix spanning all 20 samples (10 from GDS5646, comprising 5 Parkinson's and 5 healthy; 10 from GDS6063, comprising 5 influenza-exposed and 5 no-virus controls) and jointly quantile-normalized, bringing both datasets onto a shared expression-value distribution prior to any cross-condition comparison. This corrects a limitation of normalizing each dataset separately, which leaves study- and tissue-driven scale differences intact and can inflate apparent between-condition differences.

**Technical comparability and the batch–condition confound.** Formal empirical-Bayes batch-correction methods such as ComBat were deliberately **not** applied as part of the primary inferential analysis. The reason is structural rather than computational: in this study design, dataset of origin (GDS5646 versus GDS6063) and biological condition (Parkinson's versus influenza) are **perfectly confounded** — every Parkinson's sample comes from GDS5646 and every influenza sample from GDS6063. Under perfect confounding, a batch-correction model cannot distinguish variation attributable to "batch" from variation attributable to "condition": fitted without a condition covariate it would remove the very between-condition signal being tested, and fitted with condition as a covariate it is non-identifiable. Joint quantile normalization was therefore used to align marginal expression distributions, while the differential-expression results are explicitly framed below as candidate findings rather than batch-corrected, within-platform estimates. A principal-component analysis of the 20 samples before and after normalization is reported (Figure 2) as a quality-control check on residual study-level separation; because batch and condition are confounded, any dataset separation visible in this PCA cannot be attributed uniquely to biology versus technical study-of-origin effects.

---

## C. REPLACE the "Data and Code Availability" statement

> REPLACE the existing "Data and Code Availability" paragraph with:

**Data and Code Availability.** All analysis code, including data preprocessing, gene-set comparison, cross-dataset normalization, differential-expression analysis, correlation analysis, enrichment analysis, pathway mapping, and figure generation, is available in the public GitHub repository at https://github.com/nuranaverdiyeva/comparative-transcriptomic-analysis . Specifically, the end-to-end analysis is provided as a Jupyter notebook, `transcriptomic_analysis.ipynb` (https://github.com/nuranaverdiyeva/comparative-transcriptomic-analysis/blob/main/transcriptomic_analysis.ipynb), with an equivalent standalone Python script, `transcriptomic_analysis.py` (https://github.com/nuranaverdiyeva/comparative-transcriptomic-analysis/blob/main/transcriptomic_analysis.py); the figure-generation script used for the revised figures is `make_figures.py` (https://github.com/nuranaverdiyeva/comparative-transcriptomic-analysis/blob/main/make_figures.py), with all generated figures stored under the `figures/` directory (https://github.com/nuranaverdiyeva/comparative-transcriptomic-analysis/tree/main/figures). A citable, archived version of this repository is available on Zenodo (https://doi.org/10.5281/zenodo.21851016). Raw gene-expression data are publicly available from NCBI Gene Expression Omnibus under accessions GDS5646 and GDS6063 (reference series GSE54536 and GSE68849).

---

## D. REPLACE the "Functional enrichment of influenza-specific genes" results paragraph

> REPLACE the paragraph that currently begins "GO enrichment analysis of the genes unique to the influenza dataset returned more specific terms..." with:

GO enrichment analysis of the genes unique to the influenza dataset returned more specific terms than the shared-gene analysis, including several related to systemic and receptor-mediated signaling (Table 3). The appearance of *system process*, *sensory perception*, and *G protein-coupled receptor signaling* terms among genes unique to a peripheral viral-infection dataset was not anticipated and requires careful interpretation. Two non-biological explanations are plausible and are weighed against the biological reading discussed below. First, **annotation and feature-retention bias**: although both DataSets use the same platform (GPL10558), the influenza-unique set is defined by subtraction of the Parkinson's gene list, and which probes map to valid gene symbols and pass each DataSet's quality filtering differs between the two curated DataSets; sensory-perception and GPCR terms are over-represented among multi-gene receptor families, so their enrichment can arise from which probe-to-symbol annotations are retained in each DataSet rather than from coordinated biology. Second, **tissue composition**: the influenza DataSet is derived from plasmacytoid dendritic cells, a cell type with a distinctive receptor and signaling repertoire, so receptor-signaling enrichment may reflect the baseline biology of the cell type rather than an influenza-specific response. For these reasons the influenza-unique GO terms are reported descriptively and are not advanced as evidence of a shared neuroimmune mechanism; they are flagged as a finding that would need to be re-tested against a tissue- and annotation-matched background.

---

## E. REPLACE the second Discussion paragraph (the one beginning "The enrichment results unique to the influenza dataset...")

> REPLACE that paragraph with the following, which explicitly separates candidate biological signal from technical/annotation artifacts as requested by the reviewer.

The enrichment results unique to the influenza dataset, together with the differential-expression findings, are the most notable observations of this study, and we separate them here into results that are more plausibly biological and results that are more plausibly technical or annotation-driven. **More plausibly biological (hypothesis-generating only):** the differential-expression analysis independently identified several immune-signaling genes — RNF125, CD52, HCST, and ITPKB — among the most significantly different between conditions. This convergence across two independent analytic approaches (gene-set membership and per-gene differential expression) is consistent with, though not sufficient to confirm, a biological reading grounded in the innate antiviral immune literature. Viral pathogens are recognized by pattern-recognition receptors, including RIG-I-like receptors, triggering downstream NF-κB and interferon signaling essential for antiviral defense but capable of chronic inflammatory injury when persistently activated [16, 17]; RNF125 itself functions as a ubiquitin ligase that degrades RIG-I and its downstream partner MAVS, directly regulating the strength of this response [15]. Several Parkinson's-relevant genes, including LRRK2, PINK1, and PRKN, intersect with immune and mitochondrial antiviral signaling pathways, and alpha-synuclein has a proposed role in innate host defense, with viral challenge shown experimentally to increase its expression and aggregation [5, 18, 19], and a neurotropic influenza strain has been shown to enter the central nervous system and trigger neuroinflammatory and dopaminergic changes in an animal model [8]. Under this reading, the immune-signaling genes among the top differentially expressed hits could reflect systemic or neuroimmune signaling activated during infection.

**More plausibly technical or annotation-driven:** (i) the near-complete containment of the Parkinson's gene set within the influenza gene set, which most plausibly reflects DataSet-specific differences in probe-to-gene-symbol annotation and filtering rather than biology; (ii) the broad, tissue-nonspecific GO terms dominating the shared-gene enrichment (e.g., *regulation of biological process*, *signaling*), which are expected whenever the analyzed set is a large fraction of the measured transcriptome; and (iii) the influenza-unique enrichment for *sensory perception* and *G protein-coupled receptor signaling*, which can plausibly arise from the receptor-family composition of the probe annotations and from the plasmacytoid-dendritic-cell origin of the influenza DataSet, as discussed in the Results. Distinguishing a genuine neuroimmune signal from these technical and annotation effects is left as an open, testable question for follow-up work, ideally using tissue- and annotation-matched datasets with matched healthy controls and direct measurement of the specific candidate genes identified here.

---

## F. REPLACE the "Limitations" paragraph

> REPLACE the existing limitations paragraph with:

This study has several limitations. First, **batch and condition are perfectly confounded**: all Parkinson's samples come from GDS5646 and all influenza samples from GDS6063, so formal batch-correction methods such as ComBat cannot disentangle study-of-origin effects from biological condition; joint quantile normalization was used to align marginal distributions, but the differential-expression results should be interpreted as suggestive rather than as equivalent to a within-platform case-control analysis. Second, **each DataSet pools case and control samples**: GDS5646 contains 5 Parkinson's patients and 5 healthy controls, and GDS6063 contains 5 influenza-exposed and 5 no-virus samples, yet the cross-dataset comparison treats each DataSet as a single condition representative; this means part of the between-dataset difference reflects case/control status within each DataSet and tissue-of-origin (blood versus plasmacytoid dendritic cells) rather than PD-versus-influenza biology alone. Third, **sample size is small (n = 10 per DataSet; effectively n = 5 per disease subgroup)**, which limits statistical power for individual-gene-level claims and precludes replication within this study; the large effect sizes observed for the top differentially expressed genes should be interpreted cautiously pending validation in independent, matched cohorts. Fourth, this analysis is exploratory and hypothesis-generating rather than confirmatory, and the mechanistic interpretation offered above is a plausible hypothesis rather than a conclusion this dataset pairing can independently support. Future work should extend this comparison to additional Parkinson's and viral-infection datasets with matched healthy controls, apply formal batch-correction and count-based differential-expression modeling in a design where batch and condition are not confounded, and directly test the immune-signaling candidate genes identified here — particularly RNF125 and CD52 — in tissue- or platform-matched validation cohorts.

---

## G. Minor corrections

1. **Table numbering.** Renumber tables in order of first appearance: Table 1 = dataset characteristics (new, §A); Table 2 = statistical tests of gene overlap (currently "Table 1"); Table 3 = top differentially expressed genes (currently "Table 4"); Table 4 = GO terms for shared genes (currently "Table 2"); Table 5 = GO terms for influenza-unique genes (currently "Table 3"). Update all in-text references accordingly. This fixes the current anomaly whereby "Table 4" appears before "Tables 2–3".
2. **"Fifteen genes" wording.** The Results currently states "The fifteen genes with the smallest adjusted p-values are summarized in Table 4" while Table 4 lists only 11 genes. Either change "fifteen" to "eleven" (matching the 11 rows shown) or expand the table to display all 15; we recommend expanding to 15 rows (CDV3, AKTIP, GLTP, HCST are the next genes by adjusted p-value) so the table matches the sentence.
3. **KEGG reproducibility.** Soften the statement that KEGG results were drawn from "an earlier successful query window" to: "KEGG pathway mappings were obtained from the KEGG REST API during an analysis session in which the service was responsive; because the public KEGG REST endpoint intermittently returns HTTP errors, these pathway mappings should be regarded as exploratory and reproducible only when the API is available." Remove the claim that a "full re-query is planned" unless you commit to running it; if you can re-run it now (the endpoint was responsive at submission time), replace the earlier-window caveat with the current reproducible result.
4. **"Different microarray platforms" wording.** Throughout Methods and Results, replace "different microarray platforms" / "cross-platform normalization" with "different studies and tissues on the same microarray platform" / "cross-dataset normalization", because both DataSets use platform GPL10558 (Table 1). The joint normalization remains justified by cross-study and cross-tissue distributional differences.
5. **Sample labels.** Wherever the text says "10 Parkinson's samples" or "10 influenza samples" (e.g., Data sources, Methods, Results), relabel as "10 samples from GDS5646 (5 Parkinson's patients and 5 healthy controls)" and "10 samples from GDS6063 (5 influenza-exposed and 5 no-virus controls)". Figure legends and the PCA/volcano axes should likewise read "GDS5646 (PD+control)" / "GDS6063 (influenza+control)" rather than bare "PD" / "Flu", so the case/control pooling is transparent.
6. **Correlation method wording.** The Methods text describes the correlation as "mean jointly-normalized expression per gene, correlated between conditions," but the reported r = 0.009 corresponds to z-scoring each condition before averaging (a per-gene standardization step). Update the Methods *Correlation analysis* paragraph to state explicitly that each condition's expression was standardized per gene (z-scored across samples within condition) before the per-gene mean was taken and correlated, so the described method matches the reported r = 0.009. (Re-running this exact procedure on the committed data reproduces r ≈ 0.003–0.009, i.e. no coordinated expression, confirming the headline result.)
7. **Code-availability links are live only after push.** Do not paste the *Data and Code Availability* paragraph (§C) into the Cifra editor until `make_figures.py` and the `figures/` directory have been committed and pushed to GitHub; ideally cite a commit hash or an updated Zenodo version so the permalinks resolve.
8. **Gene-overlap count reconciliation.** Re-running the committed pipeline (`make_figures.py` / `transcriptomic_analysis.ipynb`) from the committed CSVs reproduces the near-complete containment finding with 12,600 shared genes, 8,162 influenza-unique, and 0 Parkinson's-unique genes, versus the 12,579 / 8,183 / 20 reported in the manuscript. The discrepancy (≈21 genes) reflects minor preprocessing/version differences and does not change any conclusion. Before final submission, re-run the canonical notebook in your environment and, if it reproduces 12,579/8,183/20, use those figures; otherwise update the manuscript text and the Table 2 contingency entries to the reproduced values so that text, tables, and figures are internally consistent.

---

## H. New figures to upload (generated from `make_figures.py`)

Upload these into the Cifra figure manager and reference them in the revised text:

- **Figure 1** — `fig2a_pca_raw.png` (panel A) and `fig2b_pca_normalized.png` (panel B) — PCA of the 20 samples before (A) and after (B) joint quantile normalization, colored by DataSet. *Cite in the new "Technical comparability and the batch–condition confound" paragraph (§B).*
- **Figure 2** — `fig3_distribution_boxplots.png` — Expression-value distributions across all 20 samples before and after joint quantile normalization. *Cite alongside Figure 1.*
- **Figure 3** — `fig4_volcano.png` — Volcano plot of per-gene Welch t-test results with Benjamini–Hochberg FDR; top genes (CALM1, CRBN, DOCK11, S100A4, CD52, RNF125, ITPKB, GPD1L) labelled. *Cite in Results, "Differential expression between conditions."*
- **Figure 4** — `fig6_go_enrichment.png` — Top GO enrichment terms for shared genes (left) and influenza-unique genes (right). *Cite in Results, "Functional enrichment."*

> **Optional QC figure (not for submission yet):** `fig1_venn_overlap.png` — a Venn diagram of gene-symbol overlap. It is held back from the submission figure list because the figure reproduces as 12,600 shared / 0 Parkinson's-unique / 8,162 influenza-unique genes from the committed CSVs, whereas the manuscript reports 12,579 / 20 / 8,183. Until you re-run the canonical `transcriptomic_analysis.ipynb` in your environment and confirm which set of counts is authoritative (see §G.8), keep the Venn as internal QC only.

Captions are provided in `reviewer_response_and_changes_summary.md`.

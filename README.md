# Biological Determinants of Adverse Drug Reaction Predictability in Drug–Gene Networks

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Mihir2728/Drug-Gene-SE-NM-Manuscript/blob/main/Drug-Gene-SE-NM-Manuscript.ipynb)

Code, data and results for the manuscript

> Chachar M.S., Yang H. *Biological Determinants of Adverse Drug Reaction Predictability in Drug–Gene Networks.* School of Mathematical & Statistical Sciences, University of Galway (manuscript in preparation).

## Overview

Most computational studies of adverse drug reactions (ADRs) aim to maximise prediction accuracy. This project asks a different question: **why are some ADRs much easier to predict than others?**

A drug similarity network is built from the drug–gene interaction profiles of 536 drugs. Then 5,599 side effects are predicted with harmonic function label propagation (Zhu et al., 2003) under five-fold cross-validation. The model is deliberately simple and interpretable: a drug's predicted score for a side effect depends only on the labels of drugs with similar gene interaction profiles. The predictability of each ADR is then explained through:

- **Target coherence:** the similarity among the drugs associated with an ADR minus their similarity to other drugs, using drug–gene (target) similarity.
- **Associated target genes:** genes over-represented among the drugs associated with each ADR.
- **KEGG pathway and GO Biological Process enrichment** of those genes, tested against the genes that drugs in the dataset actually target.

## Quick start

**Google Colab:** click the *Open in Colab* badge above, then choose **Runtime → Run all**. Nothing needs to be uploaded or mounted: the notebook downloads the two data files from this repository if they are not already present.

**Locally:**

```bash
git clone https://github.com/Mihir2728/Drug-Gene-SE-NM-Manuscript.git
cd Drug-Gene-SE-NM-Manuscript
pip install -r requirements.txt jupyter
jupyter notebook Drug-Gene-SE-NM-Manuscript.ipynb      # then Run All
```

Run time is about 10 minutes. Results are written to `outputs/` and figures to `figures/`. On Colab, the last cell also zips and downloads `outputs/`, `figures/` and `genesets/`.

## Expected results

A complete run should reproduce these values. The same values are stored in `outputs/`, so you can compare your run against them.

| Result | Value | File in `outputs/` |
|---|---|---|
| AUROC, harmonic propagation / k-NN (k = 10) | 0.670 / 0.627 | `table1_auroc_auprc_fmax.csv` |
| AUPRC, harmonic propagation / k-NN | 0.210 / 0.171 | `table1_auroc_auprc_fmax.csv` |
| F<sub>max</sub>, harmonic propagation / k-NN | 0.466 / 0.459 | `table1_auroc_auprc_fmax.csv` |
| Shuffled-network control (AUROC / AUPRC) | 0.512 / 0.107 | `network_control.csv` |
| ADRs ranked (≥ 10 drugs, all five folds) | 970 | `ranked_side_effects.csv` |
| ADRs predicted better than chance (FDR < 0.05) | 742 (76%) | `ranked_side_effects.csv` |
| Spearman ρ, target coherence vs AUROC | 0.73 | `predictability_correlations.csv` |
| ADRs with a significant KEGG pathway, best vs worst 100 | 88% vs 26% | `enrichment_group_comparison.csv` |

## Repository contents

| Path | Description |
|---|---|
| `Drug-Gene-SE-NM-Manuscript.ipynb` | Complete analysis notebook (all results, tables and figures of the manuscript) |
| `drug_sideeffect_binary_matrix.csv` | Drug × side-effect matrix (1,020 drugs × 5,599 side effects, 0/1); first column `drug_name` |
| `intersection_DGIdb_mat.tsv` | Drug × gene matrix (777 drugs × 2,022 genes, 0/1), tab-separated, drug names as row names |
| `genesets/` | The Enrichr gene-set libraries used in the manuscript (KEGG 2026, GO Biological Process 2026), saved as JSON |
| `outputs/` | All result tables produced by the notebook (see below) |
| `figures/` | All figures produced by the notebook, PNG (300 dpi) and PDF |
| `requirements.txt` | Package versions used to produce the results |
| `README.md` | This file |

### Data sources

- **Drug–gene interactions:** Drug–Gene Interaction Database (DGIdb). Any reported interaction between a drug and a gene is coded as 1.
- **Drug–side effect associations:** SIDER, as processed by Zhong et al. (2024). Side effects are MedDRA preferred terms.

These are curated database associations, not evidence that a drug causes a reaction.

## Analysis pipeline

| Notebook section | What it does | Manuscript |
|---|---|---|
| 1. Setup | Configuration; the data files are read from the repository folder or downloaded from GitHub | – |
| 2. Data | Loads both matrices; any non-zero entry is set to 1 | §2.1 |
| 3. Drug similarity and network | Jaccard similarity of gene profiles. Drugs are kept if at least one other drug has similarity ≥ 0.3 (536 of 777). The threshold only selects drugs: all non-zero similarities among the retained drugs (36,726 pairs) are used as propagation weights | §2.2 |
| 4. Harmonic prediction and k-NN | Five-fold cross-validation over drugs (`KFold`, `random_state = 27`). Harmonic solution f<sub>u</sub> = (L<sub>uu</sub> + εI)<sup>−1</sup>(−L<sub>ul</sub> f<sub>l</sub>), ε = 10<sup>−6</sup>. k-NN baseline on the same similarity matrix, k ∈ {5, 10, 15, 20} | §2.3–2.4 |
| 5. Performance summary | AUROC and AUPRC per side effect and fold. F<sub>max</sub> over thresholds 0.01–1.00 on the raw scores (no normalisation), with precision averaged over drugs with at least one prediction and recall over all drugs (Radivojac et al., 2013; Torres et al., 2021) | §2.4, Table 1 |
| 6. Predictability of individual ADRs | ADRs ranked by mean AUROC (≥ 10 drugs with and without the ADR, all five folds, identical drug sets merged). Mann–Whitney U test on out-of-fold scores with Benjamini–Hochberg correction; target coherence; shuffled-network control; similarity threshold sensitivity (0.1–0.5) | §2.5, §3.1–3.3, Table 2 |
| 7. Target genes of the top and bottom ADRs | Associated genes by hypergeometric test (BH, FDR < 0.05). KEGG and GO enrichment against the 1,017 genes targeted by the analysed drugs; gene sets with fewer than five background genes are not tested | §2.6, §3.4–3.5, Tables 3–5 |
| 8. Figures | All figures, named after their number in the manuscript | Figs. 2–6, S1–S4 |
| 9. Enrichment across the whole ranking | Enrichment of every ranked ADR. Top against bottom 10, 100 and half of the ranking (Fisher's exact test), and AUROC of ADRs enriched for each gene set (Mann–Whitney U, BH) | §3.6, Tables 6–7, Table S1 |
| 10. Package versions and download | Writes `outputs/requirements.txt` and zips the results (downloaded automatically on Colab) | – |

### Main settings (configuration cell)

| Setting | Default | Meaning |
|---|---|---|
| `DATA_DIR`, `REPO_RAW` | `"."`, this repository | Where the data files are read from, and where they are downloaded from if missing |
| `SIM_THRESHOLD` | 0.3 | Minimum gene similarity to at least one other drug for a drug to be analysed |
| `N_SPLITS`, `SEED` | 5, 27 | Cross-validation folds and random seed |
| `EPS` | 1e-6 | Regularisation of the harmonic system |
| `KNN_K_VALUES` | 5, 10, 15, 20 | k values for the baseline; the best by mean AUPRC (k = 10) is reported |
| `TIE_DECIMALS` | 10 | Scores are rounded so that values differing only by floating-point error count as ties, which makes the ranking identical across computers |
| `MIN_POSITIVE_DRUGS` | 10 | Minimum number of drugs with, and without, an ADR for it to be ranked |
| `FDR_ALPHA` | 0.05 | Benjamini–Hochberg level |
| `ENRICHR_LIBRARIES` | `KEGG_2026`, `GO_Biological_Process_2026` | Exact library versions used in the manuscript. A different year is never substituted |
| `GENESET_DIR` | `genesets` | Libraries are read from here; if absent, they are downloaded from Enrichr once and saved here |

## Outputs

### Tables (`outputs/`)

| File | Content | Manuscript |
|---|---|---|
| `table1_auroc_auprc_fmax.csv`, `fmax_per_fold.csv` | AUROC, AUPRC and F<sub>max</sub> for both models; F<sub>max</sub> and its threshold per fold | Table 1 |
| `harmonic_results_per_fold.csv`, `knn_results_per_fold.csv` | Every metric for every side effect in every fold | |
| `model_comparison.csv` | Threshold-based metrics at a fixed 0.5 cut-off (kept for reference; not used in the manuscript) | |
| `network_control.csv` | Shuffled-network control | §3.1 |
| `threshold_sensitivity.csv` | Network size and performance for thresholds 0.1–0.5 | Table 2 |
| `ranked_side_effects.csv` | The 970 ranked ADRs with AUROC, FDR, similarity measures and number of associated genes | §3.2–3.3 |
| `predictability_correlations.csv`, `predictability_quartiles.csv` | Spearman correlations and quartile summaries | §3.2–3.3 |
| `top_bottom_adr_genes.csv`, `top_bottom_adr_gene_details.csv` | Top and bottom 10 ADRs with their drugs and associated genes | Tables 3–4 |
| `top10_genes.txt`, `bottom10_genes.txt`, `background_genes.txt` | Gene lists | |
| `enrichr_results.csv`, `enrichr_summary.csv` | KEGG and GO enrichment of the top and bottom ADRs | Table 5 |
| `enrichment_per_adr.csv`, `enrichment_significant_terms_all_adrs.csv` | Enrichment of every ranked ADR | §3.6 |
| `enrichment_group_comparison.csv` | Top against bottom groups at 10, 100 and half | Table 6 |
| `pathway_predictability.csv` | Association of each gene set with predictability | Table 7, Table S1 |
| `requirements.txt` | Package versions of the run | |

### Figures (`figures/`)

| File | Manuscript |
|---|---|
| `fig2_roc_pr_curves` | Fig. 2 |
| `fig3_metrics_vs_frequency` | Fig. 3 |
| `fig4_auroc_vs_target_coherence` | Fig. 4 |
| `fig5_top_bottom_adrs` | Fig. 5 |
| `fig6_enrichment_along_ranking` | Fig. 6 |
| `figS1_gene_vs_side_effect_similarity` | Supplementary Fig. S1 |
| `figS2_network_train_test` | Supplementary Fig. S2 |
| `figS3_auroc_vs_auprc` | Supplementary Fig. S3 |
| `figS4_harmonic_vs_knn` | Supplementary Fig. S4 |
| `extra_auroc_vs_auprc_per_fold` | Not in the manuscript (per-fold version of Fig. S3) |

Fig. 1 (study workflow) was drawn separately and is not produced by the notebook.

## Reproducibility

- **Deterministic splits.** A single fixed cross-validation split (`random_state = 27`) is used throughout, so all results are deterministic.
- **Score rounding.** Harmonic scores are rounded to 10 decimals before evaluation. Without this, ADRs with nearly identical AUROC can change rank between computers.
- **Frozen gene-set libraries.** The enrichment uses the copies saved in `genesets/`, so it does not depend on the current state of the Enrichr web service. Delete a file in `genesets/` only if you want to download that library again.
- **Package versions.** The versions used for the published results are listed in `requirements.txt`.

## Limitations

- Both data sources are curated databases with annotation bias. A missing annotation does not mean a missing effect.
- 241 of the 777 drugs share too few targets with any other drug to be analysed, and ADRs reported for fewer than 10 drugs are not ranked.
- Target coherence is computed from the same drug–gene similarities used for prediction, so its association with predictability is not independent of the model. It shows how strongly the organisation of drugs in target space constrains predictability; the gene and pathway analyses provide the biological interpretation.
- Only one cross-validation partition and one type of similarity are used.

## References

- Zhu X., Ghahramani Z., Lafferty J. Semi-supervised learning using Gaussian fields and harmonic functions. *Proceedings of the 20th International Conference on Machine Learning (ICML)*, 912–919 (2003).
- Zhong Y., Seoighe C., Yang H. Non-negative matrix factorization combined with kernel regression for the prediction of adverse drug reaction profiles. *Bioinformatics Advances* 4(1), vbae009 (2024). https://doi.org/10.1093/bioadv/vbae009
- Radivojac P. et al. A large-scale evaluation of computational protein function prediction. *Nature Methods* 10(3), 221–227 (2013). https://doi.org/10.1038/nmeth.2340
- Torres M., Yang H., Romero A.E., Paccanaro A. Protein function prediction for newly sequenced organisms. *Nature Machine Intelligence* 3(12), 1050–1060 (2021).

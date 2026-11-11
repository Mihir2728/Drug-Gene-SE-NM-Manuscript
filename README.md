# Biological Determinants of Adverse Drug Reaction Predictability in Drug–Gene Networks

Code and data for the manuscript

> Chachar M.S., Yang H. *Biological Determinants of Adverse Drug Reaction Predictability in Drug–Gene Networks.* School of Mathematical & Statistical Sciences, University of Galway.

## Overview

Most computational studies of adverse drug reactions (ADRs) try to maximise prediction accuracy. This project asks a different question: **why are some ADRs much easier to predict than others?**

We build a drug similarity network from the drug–gene interaction profiles of 536 drugs. We then predict 5,599 side effects with harmonic function label propagation (Zhu et al., 2003) under five-fold cross-validation. The model is deliberately simple and interpretable: a drug's predicted score for a side effect depends only on the labels of drugs with similar gene interaction profiles. The predictability of each ADR is then explained through:

- **target coherence**: how much more the drugs associated with an ADR share gene targets with one another than with other drugs;
- **associated target genes**: genes over-represented among the drugs associated with each ADR;
- **KEGG pathway and GO Biological Process enrichment** of those genes, tested against the genes that drugs in the dataset actually target.

## Main results

| | Harmonic propagation | k-NN baseline (k = 10) |
|---|---|---|
| AUROC | 0.670 | 0.627 |
| AUPRC | 0.210 | 0.171 |
| F<sub>max</sub> | 0.466 | 0.459 |

- With the drug labels of the network shuffled, performance falls to near chance (AUROC 0.512, AUPRC 0.107, F<sub>max</sub> 0.419).
- 970 ADRs were reported for at least 10 drugs and evaluable in all folds. Of these, 742 (76%) were predicted significantly better than chance (one-sided Mann–Whitney U, Benjamini–Hochberg FDR < 0.05). Mean AUROC ranged from 0.28 to 0.96.
- Predictability is strongly associated with target coherence (Spearman's ρ = 0.73).
- The most predictable ADRs are class effects with shared receptor targets, for example:
  - coronary arteriospasm (HTR1B, HTR1D)
  - serotonin syndrome (HTR1A, HTR1B, SLC6A4)
  - tardive dyskinesia (DRD2, HTR2A)

  The least predictable are non-specific reactions reported for pharmacologically unrelated drugs.
- 88% of the 100 best predicted ADRs, but only 26% of the 100 worst, have at least one significant KEGG pathway. The distinguishing pathways are mainly neurotransmitter receptor signalling.

## Repository contents

| File | Description |
|---|---|
| `Drug-Gene-SE-NM-Manuscript.ipynb` | Complete analysis notebook (all results, tables and figures of the manuscript) |
| `drug_sideeffect_binary_matrix.csv` | Drug × side-effect matrix (1,020 drugs × 5,599 side effects, 0/1); first column `drug_name` |
| `intersection_DGIdb_mat.tsv` | Drug × gene matrix (777 drugs × 2,022 genes, 0/1), tab-separated, drug names as row names |
| `README.md` | This file |

### Data sources

- **Drug–gene interactions:** Drug–Gene Interaction Database (DGIdb). Any reported interaction between a drug and a gene is coded as 1.
- **Drug–side effect associations:** SIDER, as processed by Zhong et al. (2024). Side effects are MedDRA preferred terms.

These are curated database associations, not evidence that a drug causes a reaction.

## Requirements

Python 3.9 or later with:

```
numpy  pandas  scipy  scikit-learn  networkx  matplotlib  gseapy
```

`gseapy` is used only to download the Enrichr gene-set libraries; the notebook installs it automatically if it is missing. The enrichment step needs internet access. Everything else runs offline.

## How to run

### Google Colab (as used for the manuscript)

1. Upload `drug_sideeffect_binary_matrix.csv` and `intersection_DGIdb_mat.tsv` to `MyDrive` in Google Drive.
2. Open the notebook in Colab and choose **Runtime → Run all**. Allow access to Google Drive when asked.
3. The last cell zips `outputs/` and `figures/` and downloads them.

### Locally

```bash
pip install numpy pandas scipy scikit-learn networkx matplotlib gseapy jupyter
jupyter notebook Drug-Gene-SE-NM-Manuscript.ipynb
```

Set `DATA_DIR` in the configuration cell to the folder that contains the two data files. The Google Drive cell is skipped automatically outside Colab.

Run time is about 7–10 minutes plus a few minutes for the enrichment step.

## Analysis pipeline

| Notebook section | What it does | Manuscript |
|---|---|---|
| 2. Data | Loads both matrices; any non-zero entry is set to 1 | §2.1 |
| 3. Drug similarity and network | Jaccard similarity of gene profiles. Drugs are kept if at least one other drug has similarity ≥ 0.3 (536 of 777). The threshold only selects drugs: all non-zero similarities among the retained drugs (36,726 pairs) are used as propagation weights | §2.2 |
| 4. Harmonic prediction and k-NN | Five-fold cross-validation over drugs (`KFold`, `random_state = 27`). Harmonic solution f<sub>u</sub> = (L<sub>uu</sub> + εI)<sup>−1</sup>(−L<sub>ul</sub> f<sub>l</sub>), ε = 10<sup>−6</sup>. k-NN baseline on the same similarity matrix, k ∈ {5, 10, 15, 20} | §2.3–2.4 |
| 5. Performance summary | AUROC and AUPRC per side effect and fold. F<sub>max</sub> over thresholds 0.01–1.00 on the raw scores (no normalisation), with precision averaged over drugs with at least one prediction and recall over all drugs (Radivojac et al., 2013; Torres et al., 2021) | §2.4, Table 1 |
| 6. Predictability of individual ADRs | ADRs ranked by mean AUROC (≥ 10 drugs with and without the ADR, all five folds, identical drug sets merged). Mann–Whitney U test on out-of-fold scores with BH correction; target coherence; shuffled-network control; similarity threshold sensitivity (0.1–0.5) | §2.5, §3.1–3.3, Table 2 |
| 7. Target genes of the top and bottom ADRs | Associated genes by hypergeometric test (BH, FDR < 0.05). KEGG_2026 and GO_Biological_Process_2026 enrichment against the 1,017 genes targeted by the analysed drugs; gene sets with fewer than five background genes are not tested | §2.6, §3.4–3.5, Tables 3–5 |
| 8. Figures | All figures as PNG (300 dpi) and PDF | Figs. 2–6, S1–S4 |
| 9. Enrichment across the whole ranking | Enrichment of every ranked ADR. Top against bottom 10, 100 and half of the ranking (Fisher's exact test), and AUROC of ADRs enriched for each gene set (Mann–Whitney U, BH) | §3.6, Tables 6–7, Table S1 |

### Main settings (configuration cell)

| Setting | Default | Meaning |
|---|---|---|
| `SIM_THRESHOLD` | 0.3 | Minimum gene similarity to at least one other drug for a drug to be analysed |
| `N_SPLITS`, `SEED` | 5, 27 | Cross-validation folds and random seed |
| `EPS` | 1e-6 | Regularisation of the harmonic system |
| `KNN_K_VALUES` | 5, 10, 15, 20 | k values for the baseline; the best by mean AUPRC is reported |
| `TIE_DECIMALS` | 10 | Scores are rounded so that values differing only by floating-point error count as ties (makes the ranking identical across computers) |
| `MIN_POSITIVE_DRUGS` | 10 | Minimum number of drugs with, and without, an ADR for it to be ranked |
| `FDR_ALPHA` | 0.05 | Benjamini–Hochberg level |
| `ENRICHR_LIBRARIES` | KEGG, GO_Biological_Process | The newest Enrichr version of each is used (KEGG_2026, GO_Biological_Process_2026 for the manuscript) |

## Outputs

### Tables (`outputs/`)

| File | Content | Manuscript |
|---|---|---|
| `table1_auroc_auprc_fmax.csv`, `fmax_per_fold.csv` | AUROC, AUPRC and F<sub>max</sub> for both models; F<sub>max</sub> and its threshold per fold | Table 1 |
| `harmonic_results_per_fold.csv`, `knn_results_per_fold.csv` | Every metric for every side effect in every fold | |
| `model_comparison.csv` | All threshold-based metrics at a fixed 0.5 cut-off (kept for reference; not used in the manuscript) | |
| `network_control.csv` | Shuffled-network control | §3.1 |
| `threshold_sensitivity.csv` | Network size and performance for thresholds 0.1–0.5 | Table 2 |
| `ranked_side_effects.csv` | The 970 ranked ADRs with AUROC, FDR, similarity measures and number of associated genes | §3.2–3.3 |
| `predictability_correlations.csv`, `predictability_quartiles.csv` | Spearman correlations and quartile summaries | §3.2–3.3 |
| `top_bottom_adr_genes.csv`, `top_bottom_adr_gene_details.csv` | Top and bottom 10 ADRs with drugs and associated genes | Tables 3–4 |
| `top10_genes.txt`, `bottom10_genes.txt`, `background_genes.txt` | Gene lists | |
| `enrichr_results.csv`, `enrichr_summary.csv` | KEGG and GO enrichment of the top and bottom ADRs | Table 5 |
| `enrichment_per_adr.csv`, `enrichment_significant_terms_all_adrs.csv` | Enrichment of every ranked ADR | §3.6 |
| `enrichment_group_comparison.csv` | Top against bottom groups at 10, 100 and half | Table 6 |
| `pathway_predictability.csv` | Association of each gene set with predictability | Table 7, Table S1 |

### Figures (`figures/`)

| File | Manuscript |
|---|---|
| `figD_roc_pr_curves` | Fig. 2 |
| `figE_metrics_vs_frequency` | Fig. 3 |
| `figG_auroc_vs_target_coherence` | Fig. 4 |
| `figH_top_bottom_adrs` | Fig. 5 |
| `figI_enrichment_along_ranking` | Fig. 6 |
| `figA_gene_vs_side_effect_similarity` | Fig. S1 |
| `figB_network_train_test` | Fig. S2 |
| `figC_auroc_vs_auprc` | Fig. S3 |
| `figF_harmonic_vs_knn` | Fig. S4 |
| `figS1_auroc_vs_auprc_per_fold` | not used |

Fig. 1 (study workflow) was drawn separately and is not produced by the notebook.

## Reproducibility

- A single fixed cross-validation split (`random_state = 27`) is used throughout, so all results are deterministic.
- Harmonic scores are rounded to 10 decimals before evaluation. Without this, ADRs with nearly identical AUROC can change rank between computers.
- The Enrichr libraries are updated over time. Results in the manuscript use KEGG_2026 and GO_Biological_Process_2026; a newer library version may change the enrichment results slightly. To pin a version, give the exact library names in `ENRICHR_LIBRARIES`.

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

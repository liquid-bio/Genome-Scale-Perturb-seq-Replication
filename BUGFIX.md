# Bug log: original notebooks vs this version

The original analysis (six notebooks in `code/`, removed in this version) was rewritten into
`src/gwps.py` and seven notebooks. This file lists the problems found, how each was fixed, and how the
result changed. "Before" values are taken from the saved outputs of the original notebooks or from the
original report; "after" values come from the notebooks in `notebooks/`; "paper" values are from
Replogle et al. 2022 (Cell 185:2559).

## Summary of changed conclusions

| Analysis | Before | After | Paper |
|---|---|---|---|
| Strong perturbations, K562 genome-wide | 1,966 in one notebook, 1,948 in the others | 1,962 everywhere | 1,973 |
| CORUM complex pairs vs other pairs, median r (Fig. 2B) | 0.31 vs 0.06 | 0.60 vs 0.10 | 0.61 vs 0.10 |
| HDBSCAN perturbation clusters (Fig. 2D) | 61 with `min_cluster_size=3` (paper value 4 "did not work") | 63 with the paper's parameters | 63 (Methods) / 64 (text) |
| Perturbations compared across three screens (Fig. S4) | 1,903 / 1,945 pairs, selected by knockdown | 1,205, selected by >10 DEGs | 1,206 |
| Agreement of perturbation-perturbation correlations, K562 day 6 vs day 8 / K562 vs RPE1 (Fig. S4A/C) | r = 0.68 / 0.29 | r = 0.83 / 0.33 | 0.82 / 0.37 (cophenetic) |
| Perturbations uncorrelated between K562 and RPE1 (Fig. S4E) | "no clear function" | enriched for ribosome biogenesis (RNase MRP/P, UTP genes) | – |
| Gene expression programs (Fig. 4B) | 32 programs from all 8,248 genes; 40 perturbation clusters | 39 programs from 2,273 selected genes; 63 perturbation clusters | 38 programs, 64 clusters |
| Weak perturbations, KEGG enrichment (Fig. S2F) | report: none (code without background gave 10 terms) | 3 terms, driven by KRAB zinc-finger genes; no core pathway | – |
| Mitochondrial perturbations (Fig. 6) | 106 (K562) / 122 (RPE1, plotted as "n = 140"), chosen by name keywords; 0 mtDNA genes found | 227 (K562) / 132 (RPE1) by the paper's criteria; 13 mtDNA genes | 268 / 140 |
| Exosome knockdown raises histone and lncRNA/antisense transcripts (extension) | claimed from a scatter plot, no statistics, axis labels swapped | confirmed: exosome median z 6.7 (histone) and 27.0 (lncRNA) vs ~0, p = 7e-9 and 5e-9; EXOSC1 and EXOSC10 have weak effects | – |

Conclusions that did not change: strong perturbations are enriched for ribosome, RNA transport and
spliceosome terms; Integrator subunits form three modules with C7orf26 in the INTS10/13/14 module; ISR, UPR,
erythroid and myeloid program scores separate the expected perturbations; TMEM242 knockdown resembles ATP
synthase, not complex I, on mtDNA-encoded genes.

## Repository-level problems

| Problem | Fix |
|---|---|
| Windows paths (`os.getcwd()+'\\..\\raw data\\'`); one notebook also read `data/...` on Linux; saved outputs contained personal paths | `src/gwps.py` resolves `data/` relative to the repository (or `$GWPS_DATA`) |
| `environment.yml` was a UTF-16 Windows export with Windows builds, `prefix: D:\...`, torch pinned twice; not installable on Linux | Minimal conda-forge/bioconda `environment.yml`, tested by building it |
| README described `requirements.txt`, `main_replication.py`, `scripts/`, `results/` that did not exist | README rewritten for the actual layout |
| Three kernels, cells run out of order, many outputs stale (e.g. a printed column name that did not match the code) | Every notebook runs top to bottom with `jupyter nbconvert --execute` |
| No figure was saved (`savefig` commented out, while one cell printed "saved") | All figures written to `figures/`, intermediate tables to `results/` |
| The same loading, gene-name parsing, filtering and enrichment code was copied 3-7 times with differences | One shared module, `src/gwps.py` |
| Comments, prints and plot titles in Chinese; debug cells, empty cells, commented-out code, `!pip install` inside a notebook | English throughout; removed |
| Paper PDF committed to the repository | Removed (copyright); cited instead |

## 01. QC and strong / weak perturbations (Fig. S2)

| # | Problem | Fix | Effect |
|---|---|---|---|
| 1.1 | Strong perturbations lacked the paper's "at least 25 cells" criterion in Section 1-4; other notebooks used `num_cells_unfiltered` instead of cells passing QC | `gwps.classify`: >=50 DEGs, >=25 cells passing QC (`num_cells_filtered`), >=30% knockdown | 1,966 / 1,948 → 1,962 (paper 1,973) |
| 1.2 | Perturbations whose target was not detected (`pct_expr` NaN) were excluded; the paper keeps them ("knocked down by at least 30% or was not detected") | Included | +16 strong perturbations |
| 1.3 | Every perturbation that was neither strong nor weak was labelled `off_target`, including on-target knockdowns with 5-49 DEGs | Labels `strong`, `weak`, `intermediate`, `excluded`, `control` | Labels now mean what they say |
| 1.4 | `on_target_ess` was assigned and then overwritten with the non-targeting controls; the printed "on-target" count (109) was the number of controls | Removed | – |
| 1.5 | Enrichr was run against its default whole-genome background | Background = all genes targeted in the screen (`gwps.enrich` requires it) | Strong: unchanged (Ribosome top). Weak: 10 → 3 terms, all explained by KRAB zinc-finger genes |
| 1.6 | Fig. S2J dropped `pct_expr == -1`, i.e. the perturbations with 100% knockdown | Kept | – |

## 02. Perturbation similarity, CORUM, clustering, MDE (Fig. 2)

| # | Problem | Fix | Effect |
|---|---|---|---|
| 2.1 | Genes for correlation: top 2,319 by coefficient of variation among genes with `clean_mean > 0.25`; the paper uses the union of the top DEGs of each strong perturbation and the top 30% most variable genes with mean > 0.25 UMI | `gwps.select_genes` follows the paper; per-gene Anderson-Darling statistics are not public, so the 10 genes with the largest \|z\| per perturbation are used for the first set | 2,273 genes (paper 2,319) |
| 2.2 | The target gene was not set to 0 before computing correlations (the paper masks it) | `gwps.profiles` masks the target | – |
| 2.3 | CORUM: complexes kept if `n/N > 0.66` instead of `>= 0.66`; pairs built from gene names that collapse duplicate perturbations; the "all pairs" distribution also contained the complex pairs | `>= 0.66`; pairs of perturbations; complex pairs excluded from the background distribution; Mann-Whitney test added | median r 0.31 vs 0.06 → 0.60 vs 0.10 (paper 0.61 vs 0.10). Note: the Enrichr CORUM library has fewer complexes than CORUM 3.0 (174 kept vs 327 in the paper) |
| 2.4 | HDBSCAN `min_cluster_size` changed from the paper's 4 to 3 because 4 "only gave a few clusters"; distance `sqrt(2(1-r))` | Paper parameters on `1 - r` | 61 → 63 clusters (paper 63-64). The problem was the gene selection (2.1-2.2), not the parameter |
| 2.5 | MDE initialised randomly without a seed; not reproducible | Spectral initialisation as in the paper, fixed seed (`gwps.mde`) | Reproducible embedding |
| 2.6 | Cluster colours from `tab20` repeated across 62 clusters; noise not distinguished | Noise in grey; large clusters labelled with the best-matching CORUM complex | – |

## 03. Cross-screen comparison (Fig. S4)

| # | Problem | Fix | Effect |
|---|---|---|---|
| 3.1 | Perturbations selected by knockdown (`pct_expr <= -0.3`), not by the paper's ">10 DEGs in all three screens" | Paper criterion | 1,903 / 1,945 → 1,205 (paper 1,206) |
| 3.2 | Highly variable genes chosen separately in each screen, so the two correlation matrices were computed on different genes | One gene set: well-expressed and highly variable in all three screens | r 0.68 / 0.29 → 0.83 / 0.33 (paper 0.82 / 0.37) |
| 3.3 | Fig. S4B used the intersection of the top 4,319 CV genes of each screen (a different number from 3.2) | Genes > 0.5 UMI/cell in both screens, as in the legend | Medians 0.60 and 0.33 (paper 0.50 and 0.23) |
| 3.4 | Enrichment without background; the variable `top20_well` held 10 or 8 rows and was reused for the uncorrelated set | Background = perturbations in the comparison | Agreeing perturbations: spliceosome, RNA transport (unchanged). Disagreeing (r < 0): ribosome biogenesis, previously "no clear function" |
| 3.5 | Fig. S4F: NOPCHAP1 is called C12orf45 in the data and was silently dropped (`isin`), as were three genes absent from the essential screen | Former symbol used; missing genes printed | 24 perturbations; C1orf131, SPOUT1, UTP11 not in the screen |

## 04. Integrator (Fig. 3B)

| # | Problem | Fix | Effect |
|---|---|---|---|
| 4.1 | "Highly variable genes" recomputed as `std / mean` on z-normalized data, whose mean is close to 0; the top 2,319 were genes with mean near 0, not variable genes. Computed before removing infinite values | Genes from `gwps.select_genes` for each screen | – |
| 4.2 | RPE1 version dropped the cell-number filter | Same criteria for both screens | – |
| 4.3 | Conclusion based on visual inspection | Mean within-module vs between-module r reported | K562: INTS10/13/14 + C7orf26 within 0.74 vs 0.10 with other subunits; same pattern in RPE1 (INTS7 and INTS11 absent from RPE1) |

## 05. Gene expression programs and Fig. 4B-D

| # | Problem | Fix | Effect |
|---|---|---|---|
| 5.1 | Comment said "top 2,319 genes", code kept all 8,248 (`n_genes = min(8248, len(cv))`) | Genes from notebook 02 | 32 → 39 programs (paper 38) |
| 5.2 | Perturbation clusters recomputed here (40 clusters) instead of reusing the Fig. 2 clusters | Clusters from notebook 02 | 63 clusters × 39 programs |
| 5.3 | MDE with random initialisation, no seed; zero-variance genes produced duplicate rows | Spectral initialisation, seed, standardised rows | – |
| 5.4 | Cell 17 (t-SNE twice, not in the paper) plotted labels from a previous run | Removed | – |
| 5.5 | ISR/UPR/erythroid/myeloid scores computed only for strong perturbations; each gene z-scored across perturbations and then the scores z-scored again | Paper: mean of program genes for all perturbations with >= 25 cells, then z-normalised once | 1,948 → 10,511 perturbations scored |
| 5.6 | Gene lists hard-coded without source | Signatures = the data-driven program that best matches marker genes; overlap with the original lists reported (UPR 26/28, erythroid 27/29, myeloid 39/42, ISR 7/17 genes shared) | The original lists were close to program memberships; ISR differs most |
| 5.7 | Labels only for points above an arbitrary threshold (> 6); `texts` undefined (NameError) when no point passed | Label the 15 most outlying perturbations, as in the paper | Top ISR: AARS, EIF2B3, EIF2S1, CARS; top UPR: HSPA5, DAD1, SRP68/72; top myeloid: GATA1, LDB1 |

## 06. Mitochondria (Fig. 6, S8, 7D)

| # | Problem | Fix | Effect |
|---|---|---|---|
| 6.1 | mtDNA genes searched with `startswith('MT-')` on Ensembl IDs: 0 found, so the "nuclear" response (Fig. 6A) included the mtDNA genes | Genes matched on `gene_name` | 13 mtDNA genes; nuclear = 1,715 genes > 1 UMI/cell, as in the paper |
| 6.2 | Perturbations chosen by name keywords (`COX`, `NDU`, ...) and, in RPE1, the 122 candidates plotted under the title "n = 140" | Paper criteria: mitochondrial genes with >= 50 (K562) or >= 20 (RPE1) DEGs, >= 30 cells, >= 60% knockdown | 227 (K562) and 132 (RPE1); paper 268 and 140. MitoCarta was not available; GO terms used instead |
| 6.3 | K562 essential screen used; the paper's numbers (1,715 nuclear genes) match the genome-wide screen | Genome-wide screen | – |
| 6.4 | Nuclear features = top 500 / 1,500 most variable genes; Ward linkage on correlation-matrix rows; mtDNA genes z-scored across perturbations before correlation | Paper features; Pearson correlation of profiles; average linkage | Same-complex correlation is higher on mtDNA genes than on nuclear genes for complexes III, IV and ATP synthase in both cell lines |
| 6.5 | TMEM242 comparison was visual | Mean r reported | TMEM242 vs ATP synthase r = 0.90 (K562) / 0.72 (RPE1), vs complex I 0.46 / 0.21 |

## 07. Extension: RNA exosome, histone and lncRNA transcripts

| # | Problem | Fix | Effect |
|---|---|---|---|
| 7.1 | x and y axis labels swapped | Fixed | – |
| 7.2 | Histone list contained two non-histone genes (LINC00958, HELLPAR); lncRNA list without source | Histone signature = all HIST* genes expressed > 0.25 UMI (15 genes); lncRNA signature = data-driven program (18 of its 20 genes are in the original list) | – |
| 7.3 | Only strong perturbations scored; double z-scoring | As in 5.5 | – |
| 7.4 | Claim based on a scatter plot only | One-sided Mann-Whitney test of the 11 exosome subunits vs all other perturbations | Histone median z 6.7 vs 0.0 (p = 7e-9); lncRNA 27.0 vs -0.05 (p = 5e-9); exosome subunits rank 2-9 of 10,511 for lncRNA |
| 7.5 | KEGG on high scorers without background or adjusted-p filter | Background = all scored perturbations; adj. p < 0.05 | Top terms: RNA degradation, basal transcription factors, mRNA surveillance |

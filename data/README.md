# Data

Not tracked. Run `scripts/download_data.sh` from the repository root to download:

| File | Source | Used in |
|---|---|---|
| `K562_gwps_normalized_bulk_01.h5ad` | Replogle et al. 2022, [figshare](https://plus.figshare.com/articles/dataset/20029387) | 01, 02, 03, 04, 05, 06, 07 |
| `K562_essential_normalized_bulk_01.h5ad` | same | 01, 03 |
| `rpe1_normalized_bulk_01.h5ad` | same | 03, 04, 06 |
| `CORUM_enrichr.gmt` | Enrichr library `CORUM` | 02 |
| `GO_CC_2023.gmt` | Enrichr library `GO_Cellular_Component_2023` | 06 |

Each pseudobulk file holds, for every perturbation, the mean expression of each gene z-normalized to the
non-targeting control cells, plus per-perturbation QC (`anderson_darling_counts`, `num_cells_filtered`,
`pct_expr`, ...) and per-gene statistics (`mean` = mean UMI per cell, `std`, ...).

The paper used CORUM 3.0 (2018) and MitoCarta3.0. Neither could be downloaded from the compute environment
used here, so the Enrichr CORUM library and the GO cellular-component terms starting with "Mitochondri" are
used instead (see `BUGFIX.md`).

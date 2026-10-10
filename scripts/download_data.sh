#!/bin/bash
# Download the public inputs into data/ (about 0.6 GB). Run from the repository root.
set -e
mkdir -p data && cd data
# Replogle et al. 2022 pseudobulk files (figshare, linked from https://gwps.wi.mit.edu/)
curl -L -o K562_gwps_normalized_bulk_01.h5ad      https://ndownloader.figshare.com/files/35773217
curl -L -o K562_essential_normalized_bulk_01.h5ad https://ndownloader.figshare.com/files/35780870
curl -L -o rpe1_normalized_bulk_01.h5ad           https://ndownloader.figshare.com/files/35775512
# Enrichr gene-set libraries: CORUM complexes (notebook 02) and GO cellular component (notebook 06)
curl -L -o CORUM_enrichr.gmt "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=CORUM"
curl -L -o GO_CC_2023.gmt    "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2023"
md5sum -c - <<'MD5'
a3dfaa94ea8724217f5ecb1e14a5f0c8  K562_gwps_normalized_bulk_01.h5ad
30496767641cd2e660ee6ecb5baee132  K562_essential_normalized_bulk_01.h5ad
6f1e7d6a09e2f869759e3c4526b7f171  rpe1_normalized_bulk_01.h5ad
MD5

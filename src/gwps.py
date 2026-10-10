"""Shared helpers for the genome-wide Perturb-seq replication notebooks.

All notebooks are run from the `notebooks/` folder; data are read from `data/`
(override with the GWPS_DATA environment variable).
"""
import os
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("GWPS_DATA", ROOT / "data"))
FIG = ROOT / "figures"
RES = ROOT / "results"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7,
                     "figure.dpi": 100, "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False})
CORR_CMAP = "PuOr_r"     # correlation heatmaps: orange = positive, purple = negative (as in the paper)
EXPR_CMAP = "RdBu_r"     # z-normalized expression: red = up, blue = down

FILES = {
    "k562_gw": "K562_gwps_normalized_bulk_01.h5ad",
    "k562_ess": "K562_essential_normalized_bulk_01.h5ad",
    "rpe1": "rpe1_normalized_bulk_01.h5ad",
}

MT_GENES = ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3",
            "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB"]


def load(name):
    """Read a z-normalized pseudobulk file and add target symbol / Ensembl ID columns.

    obs names look like '10023_ZC3H18_P1P2_ENSG00000158545'.
    """
    a = ad.read_h5ad(DATA / FILES[name])
    parts = a.obs_names.str.split("_")
    a.obs["target"] = parts.str[1]
    a.obs["target_id"] = parts.str[-1]
    a.obs["is_control"] = a.obs["target"] == "non-targeting"
    return a


def knockdown_ok(obs, min_kd=0.3):
    """On-target knockdown of at least `min_kd`, or target not detected (pct_expr is NaN)."""
    return (obs["pct_expr"] <= -min_kd) | obs["pct_expr"].isna()


def classify(obs, min_degs=50, max_weak_degs=5, min_cells=25, min_kd=0.3):
    """Paper STAR Methods: strong = >=50 Anderson-Darling DEGs, >=25 cells, >=30% knockdown;
    weak = same but <5 DEGs. Everything else is labelled 'intermediate' (5-49 DEGs) or
    'excluded' (fails the cell-number or knockdown criterion)."""
    ok = ~obs["is_control"] & (obs["num_cells_filtered"] >= min_cells) & knockdown_ok(obs, min_kd)
    cls = pd.Series("excluded", index=obs.index)
    cls[ok] = "intermediate"
    cls[ok & (obs["anderson_darling_counts"] >= min_degs)] = "strong"
    cls[ok & (obs["anderson_darling_counts"] < max_weak_degs)] = "weak"
    cls[obs["is_control"]] = "control"
    return cls


def profiles(a, genes=None, mask_target=True):
    """Perturbation x gene DataFrame (gene symbols as columns, perturbation IDs as index).

    Non-finite values are set to 0 and, as in the paper, the target gene of each
    perturbation is set to 0 so that knockdown itself does not drive similarity.
    """
    genes = a.var_names if genes is None else pd.Index(genes)
    X = np.nan_to_num(np.array(a[:, genes].X, dtype=np.float32), nan=0, posinf=0, neginf=0)
    if mask_target:
        hit = a.obs["target_id"].map(pd.Series(np.arange(len(genes)), index=genes))
        rows = np.where(hit.notna())[0]
        X[rows, hit.dropna().astype(int).values] = 0
    return pd.DataFrame(X, index=a.obs_names, columns=a.var.loc[genes, "gene_name"].values)


def select_genes(a, strong_ids, min_mean=0.25, top_frac=0.3, n_top=10):
    """Approximation of the paper's 2,319 'highly variable genes' (STAR Methods, p. e10):
    union of (i) the top differentially expressed genes of every strong perturbation and
    (ii) genes with mean > 0.25 UMI/cell whose variance is in the top 30%.

    The public pseudobulk files do not contain per-gene Anderson-Darling statistics, so (i)
    uses the 10 genes with the largest |z| per perturbation. The union is restricted to
    mean > 0.25 UMI/cell, which gives a gene count close to the paper's 2,319.
    """
    v = a.var
    well = v["mean"] > min_mean
    hv = v.index[well & (v["std"] >= v.loc[well, "std"].quantile(1 - top_frac))]
    z = profiles(a[strong_ids]).abs().values
    top = v.index[np.unique(np.argsort(-z, axis=1)[:, :n_top])]
    return v.index[well & v.index.isin(hv.union(top))]


def corr_upper(c):
    """Upper-triangle values of a square correlation matrix."""
    c = np.asarray(c)
    return c[np.triu_indices_from(c, k=1)]


def read_gmt(path, suffix=None):
    """Read an Enrichr .gmt library into {term: [genes]} (optionally keep terms ending with `suffix`)."""
    sets = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if suffix is None or f[0].endswith(suffix):
            sets[f[0]] = [g for g in f[2:] if g]
    return sets


def enrich(genes, background, library="KEGG_2021_Human", padj=0.05, top=10):
    """Over-representation test (gseapy/Enrichr library, hypergeometric) against `background`,
    which should be all genes targeted in the screen, not the whole genome."""
    import gseapy as gp
    res = gp.enrichr(gene_list=list(genes), gene_sets=library, background=list(background),
                     organism="human", outdir=None).results
    res = res[res["Adjusted P-value"] < padj].sort_values("Adjusted P-value").head(top).copy()
    res["-log10 adj. P"] = -np.log10(res["Adjusted P-value"])
    return res


def enrich_barplot(res, title, name):
    n = len(res)
    fig, ax = plt.subplots(figsize=(5, 0.18 * n + 0.8))   # fixed width, height per term: same bar thickness
    ax.barh(np.arange(n), res["-log10 adj. P"], color="#00327D", height=0.7)
    ax.set_yticks(np.arange(n), res["Term"])
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlabel("-log10 adjusted P")
    ax.set_title(title, loc="left")
    savefig(fig, name)


def triangle_heatmap(c, labels, colors, title, name, vmax=1):
    """Lower-triangle correlation heatmap ordered by average-linkage clustering, with a side
    colour bar for `labels` (dict label -> colour in `colors`), as in the paper's Fig. 6A/C."""
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform
    d = squareform(np.clip(1 - c, 0, 2), checks=False)
    order = leaves_list(linkage(d, method="average"))
    c = c[np.ix_(order, order)]
    fig = plt.figure(figsize=(5, 4.6))
    ax_bar = fig.add_axes([0.06, 0.08, 0.03, 0.82])
    ax = fig.add_axes([0.10, 0.08, 0.72, 0.82])
    ax_cb = fig.add_axes([0.86, 0.55, 0.03, 0.3])
    masked = np.ma.masked_array(c, mask=np.triu(np.ones_like(c, dtype=bool), k=1))
    im = ax.imshow(masked, cmap=CORR_CMAP, vmin=-vmax, vmax=vmax, interpolation="nearest", aspect="auto")
    ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]
    ax.set_title(title, loc="left")
    lab = np.asarray(labels)[order]
    cats = list(colors)
    ax_bar.imshow(np.array([[cats.index(x)] for x in lab]), aspect="auto", interpolation="nearest",
                  cmap=ListedColormap(list(colors.values())),
                  vmin=0, vmax=len(cats) - 1)
    ax_bar.set_xticks([]); ax_bar.set_yticks([])
    fig.colorbar(im, cax=ax_cb, label="Pearson r")
    handles = [plt.Rectangle((0, 0), 1, 1, color=v) for k, v in colors.items() if k in set(lab)]
    fig.legend(handles, [k for k in colors if k in set(lab)], loc="lower right", bbox_to_anchor=(0.99, 0.1), frameon=False)
    savefig(fig, name)


def savefig(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.show()


def mde(data, dim, seed=0):
    """pyMDE embedding with the paper's settings: preserve_neighbors (n_neighbors=7,
    repulsive_fraction=5), initialized with a spectral embedding (n_neighbors=7, arpack)."""
    import pymde
    import torch
    from sklearn.manifold import SpectralEmbedding
    pymde.seed(seed)
    data = np.asarray(data, dtype=np.float64)
    init = SpectralEmbedding(n_components=dim, affinity="nearest_neighbors", n_neighbors=7,
                             eigen_solver="arpack", random_state=seed).fit_transform(data)
    problem = pymde.preserve_neighbors(torch.tensor(data, dtype=torch.float32), embedding_dim=dim, n_neighbors=7, init="random",
                                       repulsive_fraction=5)
    X0 = torch.tensor((init - init.mean(0)) / init.std(0), dtype=torch.float32)  # centred, unit variance
    return problem.embed(X=X0).numpy()


def standardize_rows(M):
    """z-normalize each row (profile), so that Euclidean distance reflects correlation."""
    M = np.asarray(M, dtype=np.float64)
    sd = M.std(1, keepdims=True)
    return (M - M.mean(1, keepdims=True)) / np.where(sd == 0, 1, sd)


def cluster_names(targets, labels, complexes, min_shared=3, max_size=100):
    """Name each perturbation cluster after the best-matching CORUM complex if it shares >= min_shared
    targets, otherwise after two of its members. Very large 'complexes' (> max_size subunits) are not
    used as names. Returns {cluster: name}."""
    complexes = {c: set(g) for c, g in complexes.items() if len(set(g)) <= max_size}
    names = {}
    for k in sorted(set(labels) - {-1}):
        members = set(np.asarray(targets)[np.asarray(labels) == k])
        # score = shared^2 / (|cluster| * |complex|): favours complexes that cover the cluster and are covered by it
        best = max(complexes, key=lambda c: len(members & complexes[c]) ** 2 / len(complexes[c]))
        n = len(members & complexes[best])
        names[k] = best.replace(" (human)", "")[:32] if n >= min_shared else ", ".join(sorted(members)[:2])
    return names


def complex_library():
    """Complex annotations used to name clusters: GO cellular component terms and CORUM human complexes."""
    import re
    lib = {re.sub(r" \(GO:\d+\)", "", k): v for k, v in read_gmt(DATA / "GO_CC_2023.gmt").items()}
    lib.update(read_gmt(DATA / "CORUM_enrichr.gmt", suffix="(human)"))
    return lib

from typing import Sequence
from typing import Union

import anndata as ad
import numpy as np
import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/common/cxg_mouse_pancreas_atlas/dataset.h5ad",
    "alpha": 1,
    "n_obs": 100,
    "cell_lb": 10,
    "cell_ub": 30,
    "umi_lb": 1000,
    "umi_ub": 5000,
    "simulated_data": "dataset_simulated.h5ad"
}
meta = {
    "functionality_name": "dataset_simulator",
    "resources_dir": "src/tasks/spatial_decomposition/dataset_simulator",
}
## VIASH END

CELLTYPE_MIN_CELLS = 25

# Reading input dataset
adata = ad.read_h5ad(par['input'])


def generate_synthetic_dataset(
    adata: ad.AnnData,
    alpha: Union[float, Sequence] = 1.0,
    n_obs: int = 1000,
    cell_lb: int = 10,
    cell_ub: int = 30,
    umi_lb: int = 1000,
    umi_ub: int = 5000,
) -> ad.AnnData:
    """Create cell-aggregate samples for ground-truth spatial decomposition task.

    Parameters
    ----------
    adata: AnnData
        Anndata object.
    type_column: str
        name of column in `adata.obs` where cell type labels are given
    alpha: Union[float,Sequence]
        alpha value in dirichlet distribution. If single number then all alpha_i values
        will be set to this value. Default value is 1.
    n_obs: int
        number of spatial observations to generate. Default value is 1000.
    cell_lb: int
        lower bound for number of cells at each spot. Default value is 10.
    cell_ub: int
        upper bound for number of cells at each spot. Default value is 30.
    umi_lb: int
        lower bound for number of UMIs at each spot. Default value is 10.
    umi_ub: int
        upper bound for number of UMIs at each spot. Default value is 30.

    Returns
    -------
    AnnData with:
        - `adata_merged.X`: simulated counts (aggregate of sc dataset).
        - `adata_merged.obsm["proportions_true"]`: true proportion values.
        - `adata_merged.obsm["coordinates"]`: coordinates of each spot.
        - `adata_merged.obsm["n_cells"]`: number of cells from each type at every location.

    """
    
    # remove rare celltypes
    adata = filter_celltypes(adata)

    # set random generator seed
    rng = np.random.default_rng(42)

    # get single cell expression data
    counts = adata.layers['counts']
    # get cell annotations/labels
    labels = adata.obs['cell_type'].values
    # get unique labels
    uni_labs = np.unique(labels)
    # count number of labels
    n_labs = len(uni_labs)
    # get number of genes
    n_genes = adata.shape[1]

    # create dict with indices of each label
    label_indices = dict()
    for label in uni_labs:
        label_indices[label] = np.where(labels == label)[0]

    # adjust alpha to vector if single scalar
    if not hasattr(alpha, "__len__"):
        alpha = np.ones(n_labs) * alpha
    else:
        assert len(alpha) == n_labs, "alpha must be same size as number of cell types"

    # generate probability of sampling label at each spot
    sp_props = rng.dirichlet(alpha, size=n_obs)
    # number of cells present at each spot
    n_cells = rng.integers(cell_lb, cell_ub, size=n_obs)

    # initialize spatial expression matrix
    sp_x = np.zeros((n_obs, n_genes))
    # initialize spatial proportion matrix
    sp_p = np.zeros((n_obs, n_labs))
    # initialize spatial cell number matrix
    sp_c = np.zeros(sp_p.shape)

    # generate expression vector for each spot (s)
    for s in range(n_obs):
        # number of cells from each label at s
        raw_s = rng.multinomial(n_cells[s], pvals=sp_props[s, :])
        # store number of cells from each type at s
        sp_c[s, :] = raw_s
        # compute proportion of each type at s
        prop_s = raw_s / n_cells[s]
        # store proportion of each type at s
        sp_p[s, :] = prop_s

        # initialize transcript pool at s
        pool_s = np.zeros(n_genes)

        # add molecules to transcript pool
        for lab, n in enumerate(raw_s):
            # get indices of cells from which transcripts should be added
            idx_sl = rng.choice(label_indices[uni_labs[lab]], size=n)
            # add molecules to pool
            pool_s += counts[idx_sl, :].sum(axis=0).A.flatten()

        # number of UMIs at spot s
        n_umis = rng.integers(umi_lb, umi_ub)
        # compute probability of sampling UMI from gene
        prob_pool_s = pool_s / pool_s.sum()

        # sample transcripts from pool
        sp_x[s, :] = np.random.multinomial(n=n_umis, pvals=prob_pool_s)

    obs_names = ["spatial_{}".format(x) for x in range(n_obs)]
    adata_spatial = ad.AnnData(
        sp_x,
        obs=dict(obs_names=obs_names),
        var=dict(var_names=adata.var_names),
    )

    # fake coordinates
    adata_spatial.obsm["coordinates"] = rng.random((adata_spatial.shape[0], 2))
    adata_spatial.obsm["proportions_true"] = sp_p
    adata_spatial.obs["n_cells"] = n_cells
    adata_spatial.obsm["n_cells"] = sp_c
    
    adata_merged = ad.concat(
        {"sc": adata, "sp": adata_spatial}, 
        label="modality",
        join="outer", 
        index_unique=None, 
        merge="unique", 
        uns_merge="unique"
    )
    adata_merged.X[adata_merged.X == np.inf] = adata_merged.X.max()  # remove inf
    adata_merged.layers["counts"] = adata_merged.X
    adata_merged.uns["cell_type_names"] = uni_labs
    return adata_merged


def filter_celltypes(adata, min_cells=CELLTYPE_MIN_CELLS):
    """Filter rare celltypes from an AnnData"""
    celltype_counts = adata.obs["cell_type"].value_counts() >= min_cells
    keep_cells = np.isin(adata.obs["cell_type"], celltype_counts.index[celltype_counts])
    return adata[adata.obs.index[keep_cells]].copy()


def filter_genes_cells(adata):
    """Remove empty cells and genes."""
    if "var_names_all" not in adata.uns:
        # fill in original var names before filtering
        adata.uns["var_names_all"] = adata.var.index.to_numpy()
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=2)


adata.X = adata.layers["counts"]
sc.pp.filter_genes(adata, min_counts=10)
adata_merged = generate_synthetic_dataset(adata, 
    alpha=par['alpha'], 
    n_obs=par['n_obs'], 
    cell_lb=par['cell_lb'], 
    cell_ub=par['cell_ub'], 
    umi_lb=par['umi_lb'], 
    umi_ub=par['umi_ub'] 
)
adata_merged.uns["spatial_data_summary"] = f"Dirichlet alpha={par['alpha']}"
filter_genes_cells(adata_merged)
adata_merged.X = None

# Convert non-string objects to categoricals to avoid
# TypeError: Can't implicitly convert non-string objects to strings
# In this case, the error is raised when there are NA values in .obs columns with dtype object (boolean).
# The resulting anndata object cannot be written to a file.
# This conversion is handled in later versions of anndata (0.10)
for col in adata_merged.obs:
    if adata_merged.obs[col].dtype == 'object':
        adata_merged.obs[col] = adata_merged.obs[col].astype('category')

print("Writing output to file")
adata_merged.write_h5ad(par["simulated_data"])

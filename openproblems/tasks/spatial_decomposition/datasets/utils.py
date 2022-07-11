from ..utils import merge_sc_and_sp
from typing import Sequence
from typing import Union

import anndata as ad
import numpy as np


def generate_synthetic_dataset(
    adata: ad.AnnData,
    type_column: str = "label",
    alpha: Union[float, Sequence] = 1.0,
    n_obs: int = 1000,
    cell_lb: int = 10,
    cell_ub: int = 30,
    umi_lb: int = 1000,
    umi_ub: int = 5000,
    seed: int = 42,
) -> ad.AnnData:
    """Create cell-aggregate samples for ground-truth spatial decomposition task.

    Parameters
    ----------
    adata: AnnData
        Anndata object.
    type_column: str
        name of column in `adata.obs` where cell type labels are gives
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
    seed: int
        Seed for rng.

    Returns
    -------
    AnnData with:
        - `adata_spatial.X`: simulated counts (aggregate of sc dataset).
        - `adata_spatial.uns["sc_reference"]`: original sc adata for reference.
        - `adata_spatial.obsm["proportions_true"]`: true proportion values.
        - `adata_spatial.obsm["n_cells"]`: number of cells from each type at
           every location
        - `adata_spatial.obs["proportions_true"]`:
           total number of cells at each location

    The cell type labels are stored in adata_sc.obs["label"].
    """

    # set random generator seed
    rng = np.random.default_rng(seed)

    # get single cell expression data
    X = adata.X
    # get cell annotations/labels
    labels = adata.obs[type_column].values
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
            pool_s += X[idx_sl, :].sum(axis=0)

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
    adata_spatial.obsm["spatial"] = rng.random((adata_spatial.shape[0], 2))
    adata_spatial.obsm["proportions_true"] = sp_p
    adata_spatial.obs["n_cells"] = n_cells
    adata_spatial.obsm["n_cells"] = sp_c
    adata_merged = merge_sc_and_sp(adata, adata_spatial)
    adata_merged.X[adata_merged.X == np.inf] = adata_merged.X.max()  # remove inf
    adata_merged.layers["counts"] = adata_merged.X.copy()

    return adata_merged

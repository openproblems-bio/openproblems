"""
From:
https://github.com/romain-lopez/DestVI-reproducibility/blob/master/simulations/
"""
from ...utils import merge_sc_and_sp
from numba import jit
from pathlib import Path
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph
from typing import Optional

import anndata
import numpy as np
import pandas as pd


def categorical(p, n_samples):
    size = list(p.shape[:-1])
    size.insert(0, n_samples)
    return (p.cumsum(-1) >= np.random.uniform(size=size)[..., None]).argmax(-1).T


@jit(nopython=True)
def get_mean_normal(cell_types, gamma, mean_, components_):  # pragma: no cover
    """Util for preparing the mean of the normal distribution.

    cell_types: (n_spots, n_cells)
    gamma: (n_spots, n_cells, n_latent)

    return: samples: (n_spots, n_cells, n_genes)
    """
    # extract shapes
    n_spots = gamma.shape[0]
    n_cells = gamma.shape[1]
    n_genes = components_[0].shape[1]

    mean_normal = np.zeros((n_spots, n_cells, n_genes))
    for spot in range(n_spots):
        for cell in range(n_cells):
            mean_normal[spot, cell] = mean_[cell_types[spot, cell]]
            c = components_[cell_types[spot, cell]]
            g = np.expand_dims(gamma[spot, cell], 0)
            mean_normal[spot, cell] += np.dot(g, c)[0]
    return mean_normal


def generate_synthetic_dataset(
    test: bool,
    lam_ct: float = 0.1,
    temp_ct: float = 1.0,
    lam_gam: float = 0.5,
    sf_gam: float = 15.0,
    bin_sampling: float = 1.0,
    ct_study: int = 0,
    grid_size: Optional[int] = None,  # size of spatial grid
    K_sampled: Optional[int] = None,  # cells sampled for each spot
    seed: int = 0,
):
    import torch

    np.random.seed(seed)
    torch.manual_seed(seed)

    # parameters
    K = 100
    if test:
        K_sampled = K_sampled or 3
        grid_size = grid_size or 5
    else:  # pragma: nocover
        K_sampled = K_sampled or 20
        grid_size = grid_size or 10

    script_dir = Path(__file__).resolve().parent

    grtruth_PCA = np.load(script_dir.joinpath("ground_truth_pca.npz").as_posix())
    mean_, components_ = grtruth_PCA["mean_"], grtruth_PCA["components_"]

    inv_dispersion = np.load(script_dir.joinpath("inverse_dispersion.npy").as_posix())

    C = components_.shape[0]
    D = components_.shape[1]

    locations, freq_sample, gamma = generate_spatial_information(
        C=C,
        D=D,
        grid_size=grid_size,
        lam_ct=lam_ct,
        temp_ct=temp_ct,
        lam_gam=lam_gam,
        sf_gam=sf_gam,
    )

    cell_types_sc = categorical(freq_sample, K)
    gamma_sc = gamma[:, None, :].repeat(K, axis=1)
    location_sc = locations[:, None, :].repeat(K, axis=1)
    # get means of the Gaussian using the sPCA model
    mean_normal = get_mean_normal(cell_types_sc, gamma_sc, mean_, components_)
    # convert back to count distribution and sample from Poisson
    mean_normal[mean_normal <= 0] = np.min(mean_normal[mean_normal > 0]) * 0.01
    transformed_mean = np.expm1(mean_normal)

    # dispersion was learned on the single-cell data.
    # this simulation might have different library sizes
    # we must match them so that the dispersion estimates make sense
    # (i.e., genes are as overpoissonian as in the experiments)
    inv_dispersion *= 1e2

    # Important remark: Gamma is parametrized by the rate = 1/scale!
    gamma_s = torch.distributions.Gamma(
        concentration=torch.tensor(inv_dispersion),
        rate=torch.tensor(inv_dispersion) / torch.tensor(transformed_mean),
    ).sample()
    mean_poisson = torch.clamp(gamma_s, max=1e8).cpu().numpy()
    transformed_mean = mean_poisson

    samples = np.random.poisson(lam=transformed_mean)

    X = samples[:, :K_sampled].reshape((-1, samples.shape[-1]))
    obs_names = [f"sc_{i}" for i in range(X.shape[0])]
    sc_anndata = anndata.AnnData(
        X=X,
        obs=dict(obs_names=obs_names),
    )
    sc_anndata.obs["cell_type"] = cell_types_sc[:, :K_sampled].reshape(-1, 1)
    sc_anndata.obs["label"] = sc_anndata.obs["cell_type"].astype(str).astype("category")
    sc_anndata.obs["n_counts"] = np.sum(sc_anndata.X, axis=1)
    sc_anndata.obsm["gamma"] = gamma_sc[:, :K_sampled].reshape(-1, gamma.shape[-1])
    sc_anndata.obsm["locations"] = location_sc[:, :K_sampled].reshape(-1, 2)

    # cluster the single-cell data using sklearn
    target_list = [2, 4, 8, 16]
    key_list = ["cell_type"]
    hier_labels_sc = np.zeros((sc_anndata.n_obs, len(target_list)))
    for ct in range(C):
        slice_ind = np.where(sc_anndata.obs["cell_type"] == ct)
        slice_counts = sc_anndata.X[slice_ind]
        slice_normalized = slice_counts / np.sum(slice_counts, axis=1)[:, np.newaxis]
        slice_embedding = PCA(n_components=10).fit_transform(
            np.log(1 + 1e4 * slice_normalized)
        )
        knn_graph = kneighbors_graph(slice_embedding, 30, include_self=False)
        for i, target in enumerate(target_list):
            labels = AgglomerativeClustering(
                n_clusters=target, connectivity=knn_graph
            ).fit_predict(slice_embedding)
            hier_labels_sc[slice_ind, i] = labels

    # aggregate hierarchical labels and append to anndata
    for i, target in enumerate(target_list):
        base_cell_type = sc_anndata.obs["cell_type"]
        sub_cell_type = hier_labels_sc[:, i]
        nb_sub_ct = len(np.unique(sub_cell_type))
        all_cell_type = np.array(
            [
                base_cell_type[j] * nb_sub_ct + sub_cell_type[j]
                for j in range(sc_anndata.n_obs)
            ]
        )
        key = str(target) + "th_sub-cell_type"
        sc_anndata.obs[key] = all_cell_type.astype(np.int)
        key_list.append(key)
    # dump keys as well
    sc_anndata.uns["key_clustering"] = key_list
    sc_anndata.uns["target_list"] = [1] + target_list

    # write the full data
    # sc_anndata.write(output_dir + "sc_simu.h5ad", compression="gzip")
    # remove the "last" cell type and dump to separate file
    # sc_anndata_partial = sc_anndata[sc_anndata.obs["cell_type"] != C - 1]
    # sc_anndata_partial.write(output_dir + "sc_simu_partial.h5ad", compression="gzip")

    transformed_mean_st_full = transformed_mean.mean(1)
    # here we should create a mask to remove contributions from one cell type
    transformed_mean_st_partial = np.sum(
        (cell_types_sc != C - 1)[:, :, np.newaxis] * transformed_mean, axis=1
    )
    transformed_mean_st_partial /= np.sum(cell_types_sc != C - 1, 1)[:, np.newaxis]

    if ct_study == 1:
        list_transformed = [transformed_mean_st_full, transformed_mean_st_partial]
    elif ct_study == 0:
        list_transformed = [transformed_mean_st_full]
    else:
        raise NotImplementedError
    for i, transformed_mean_st in enumerate(list_transformed):
        # Important remark: Gamma is parametrized by the rate = 1/scale!
        gamma_st = torch.distributions.Gamma(
            concentration=torch.tensor(inv_dispersion),
            rate=torch.tensor(inv_dispersion) / torch.tensor(transformed_mean_st),
        ).sample()
        mean_poisson_st = torch.clamp(gamma_st, max=1e8).cpu().numpy()
        mean_st = mean_poisson_st

        samples_st = np.random.poisson(lam=mean_st)
        samples_st = np.random.binomial(samples_st, bin_sampling)
        obs_names = [f"st_{i}" for i in range(samples_st.shape[0])]
        st_anndata = anndata.AnnData(
            X=samples_st,
            obs=dict(obs_names=obs_names),
        )
        if i == 0:
            altered_freq = freq_sample
        if i > 0:
            altered_freq = freq_sample[:, :-1]
            altered_freq /= altered_freq.sum(1)[:, np.newaxis]
        columns = np.arange(altered_freq.shape[1]).astype(np.str)
        freq_df = pd.DataFrame(
            altered_freq, index=st_anndata.obs.index, columns=columns
        )
        # st_anndata.obsm["cell_type"] = freq_df
        st_anndata.obsm["proportions_true"] = freq_df
        st_anndata.obsm["gamma"] = gamma
        st_anndata.obsm["locations"] = locations
        st_anndata.obs["n_counts"] = np.sum(st_anndata.X, axis=1)
        st_anndata.uns["key_clustering"] = key_list
        st_anndata.uns["target_list"] = [1] + target_list

    sc_anndata.layers["counts"] = sc_anndata.X.copy()
    st_anndata.layers["counts"] = st_anndata.X.copy()

    merged_anndata = merge_sc_and_sp(sc_anndata, st_anndata)

    return merged_anndata


def generate_spatial_information(
    grid_size,
    C,
    lam_ct,
    temp_ct,
    lam_gam,
    D,
    sf_gam,
):
    locations = (
        np.mgrid[-grid_size:grid_size:0.5, -grid_size:grid_size:0.5].reshape(2, -1).T
    )
    # get the kernel bandwidth for GP simulation
    dist_table = pdist(locations)
    bandwidth = np.median(dist_table)
    # sample from the multivariate GP for cell type
    K = np.exp(-squareform(dist_table) ** 2 / (lam_ct * bandwidth**2))
    N = K.shape[0]
    sample = np.random.multivariate_normal(np.zeros(N), K, size=C).T
    # get through softmax
    e_sample = np.exp(sample / temp_ct)
    freq_sample = e_sample / np.sum(e_sample, 1)[:, np.newaxis]

    # form the multivariate GP covariance for gamma
    K = sf_gam * np.exp(-squareform(dist_table) ** 2 / (lam_gam * bandwidth**2))
    # get latent variable for each cell types
    gamma = np.random.multivariate_normal(np.zeros(N), K, size=(D)).T

    return locations, freq_sample, gamma

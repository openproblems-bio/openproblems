from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph
from torch.distributions import Gamma
from utils import categorical
from utils import get_mean_normal

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

np.random.seed(0)


torch.manual_seed(0)


def generate_synthetic_dataset(
    input_file: str = "_input_data/",
    lam_ct: float = 0.1,
    temp_ct: float = 1.0,
    lam_gam: float = 0.5,
    sf_gam: float = 15.0,
    bin_sampling: float = 1.0,
    ct_study: int = 1,
):
    # parameters
    K = 100
    K_sampled = 20
    grid_size = 10

    grtruth_PCA = np.load(f"{input_file}grtruth_PCA.npz")
    mean_, components_ = grtruth_PCA["mean_"], grtruth_PCA["components_"]

    inv_dispersion = np.load(f"{input_file}inv-dispersion.npy")

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
        savefig=False,
    )

    # put together a summary of distinct cell types per spot
    cell_types_sc = categorical(freq_sample, 10)
    x = np.sort(cell_types_sc, axis=1)
    res = (x[:, 1:] != x[:, :-1]).sum(axis=1) + 1
    plt.hist(res, bins=np.linspace(0, C, 20))
    plt.xlabel("Number of cell types")
    plt.ylabel("Number of spots")
    plt.title(f"temp-ct={temp_ct}")
    plt.tight_layout()
    # plt.savefig(output_dir+"fre.png")

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

    if True:
        # Important remark: Gamma is parametrized by the rate = 1/scale!
        gamma_s = Gamma(
            concentration=torch.tensor(inv_dispersion),
            rate=torch.tensor(inv_dispersion) / torch.tensor(transformed_mean),
        ).sample()
        mean_poisson = torch.clamp(gamma_s, max=1e8).cpu().numpy()
        transformed_mean = mean_poisson

    samples = np.random.poisson(lam=transformed_mean)

    sc_anndata = anndata.AnnData(
        X=csr_matrix(samples[:, :K_sampled].reshape((-1, samples.shape[-1])))
    )
    sc_anndata.obs["cell_type"] = cell_types_sc[:, :K_sampled].reshape(-1, 1)
    sc_anndata.obs["label"] = pd.Categorical(
        sc_anndata.obs["cell_type"].values.astype(str)
    )
    sc_anndata.obs["n_counts"] = np.sum(sc_anndata.X.A, axis=1)
    sc_anndata.obsm["gamma"] = gamma_sc[:, :K_sampled].reshape(-1, gamma.shape[-1])
    sc_anndata.obsm["locations"] = location_sc[:, :K_sampled].reshape(-1, 2)

    # cluster the single-cell data using sklearn
    target_list = [2, 4, 8, 16]
    key_list = ["cell_type"]
    hier_labels_sc = np.zeros((sc_anndata.n_obs, len(target_list)))
    for ct in range(C):
        slice_ind = np.where(sc_anndata.obs["cell_type"] == ct)
        slice_counts = sc_anndata.X[slice_ind].A
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
    # file_name = ["st_simu.h5ad", "st_simu_partial.h5ad"]
    for i, transformed_mean_st in enumerate(list_transformed):
        # Important remark: Gamma is parametrized by the rate = 1/scale!
        gamma_st = Gamma(
            concentration=torch.tensor(inv_dispersion),
            rate=torch.tensor(inv_dispersion) / torch.tensor(transformed_mean_st),
        ).sample()
        mean_poisson_st = torch.clamp(gamma_st, max=1e8).cpu().numpy()
        mean_st = mean_poisson_st

        samples_st = np.random.poisson(lam=mean_st)
        samples_st = np.random.binomial(samples_st, bin_sampling)

        st_anndata = anndata.AnnData(X=csr_matrix(samples_st))
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
        st_anndata.obsm["n_counts"] = np.sum(st_anndata.X, axis=1)
        st_anndata.uns["key_clustering"] = key_list
        st_anndata.uns["target_list"] = [1] + target_list
        st_anndata.uns["sc_reference"] = sc_anndata
        # st_anndata.write(output_dir + file_name[i], compression="gzip")
        if i == 0:
            plt.figure(figsize=(5, 5))
            plt.hist(st_anndata.obsm["n_counts"], bins=100)
            plt.xlabel("Number of UMIs")
            plt.ylabel("Number of spots")
            plt.title(f"bin-sampling={bin_sampling}")
            plt.tight_layout()
            # plt.savefig(output_dir+"lib.png")
    return st_anndata


def generate_spatial_information(
    grid_size, C, lam_ct, temp_ct, lam_gam, D, sf_gam, savefig: bool
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

    # plot cell types
    if savefig:
        plt.figure(figsize=(9, 9))
        for i in range(0, C):
            plt.subplot(331 + i)
            plt.scatter(locations[:, 0], locations[:, 1], c=freq_sample[:, i])
            plt.title("cell type " + str(i))
            plt.colorbar()
        plt.tight_layout()
        plt.savefig(savefig + "cell_type_proportion.png")
        plt.clf()
        # plot gammas
        plt.figure(figsize=(5, 5))
        for i in range(D):
            plt.subplot(221 + i)
            plt.scatter(locations[:, 0], locations[:, 1], c=gamma[:, i])
            plt.title("gamma " + str(i))
            plt.colorbar()
        plt.tight_layout()
        plt.savefig(savefig + "gamma.png")
        plt.clf()
    return locations, freq_sample, gamma

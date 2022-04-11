from ....tools.normalize import log_cpm

import scanpy as sc


def preprocess_scanpy(adata):

    adata = log_cpm(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    adata.obsm["X_input"] = adata.obsm["X_pca"]

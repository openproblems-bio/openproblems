from .....tools.decorators import baseline_method
from .....tools.utils import check_version
from ...batch_integration_graph.methods.baseline import _random_embedding
from ...batch_integration_graph.methods.baseline import _randomize_features

import numpy as np
import scanpy as sc


@baseline_method(
    method_name="No Integration",
    method_summary="TODO",
)
def no_integration(adata, test=False):
    adata.obsm["X_emb"] = adata.obsm["X_uni_pca"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Integration",
    method_summary="TODO",
)
def random_integration(adata, test=False):
    adata.obsm["X_emb"] = _randomize_features(adata.obsm["X_uni_pca"])
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Integration by Celltype",
    method_summary="TODO",
)
def celltype_random_integration(adata, test=False):
    adata.obsm["X_emb"] = _randomize_features(
        adata.obsm["X_uni_pca"], partition=adata.obs["labels"]
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Embedding by Celltype",
    method_summary="TODO",
)
def celltype_random_embedding(adata, test=False):
    adata.obsm["X_emb"] = _random_embedding(partition=adata.obs["labels"])
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Integration by Batch",
    method_summary="TODO",
)
def batch_random_integration(adata, test=False):
    adata.obsm["X_emb"] = _randomize_features(
        adata.obsm["X_uni_pca"], partition=adata.obs["batch"]
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="No Integration by Batch",
    method_summary="TODO",
)
def no_integration_batch(adata, test=False):
    """Compute PCA independently on each batch

    See https://github.com/theislab/scib/issues/351
    """
    adata.obsm["X_emb"] = np.zeros((adata.shape[0], 50), dtype=float)
    for batch in adata.obs["batch"].unique():
        batch_idx = adata.obs["batch"] == batch
        n_comps = min(50, np.sum(batch_idx))
        solver = "full" if n_comps == np.sum(batch_idx) else "arpack"
        adata.obsm["X_emb"][batch_idx, :n_comps] = sc.tl.pca(
            adata[batch_idx],
            n_comps=n_comps,
            use_highly_variable=False,
            svd_solver=solver,
            copy=True,
        ).obsm["X_pca"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

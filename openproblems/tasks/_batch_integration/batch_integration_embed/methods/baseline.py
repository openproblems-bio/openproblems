from .....tools.decorators import baseline_method
from .....tools.utils import check_version
from ..._common.methods.baseline import _random_embedding

import numpy as np
import scanpy as sc


@baseline_method(
    method_name="Random Embedding by Celltype (with jitter)",
    method_summary=(
        "Cells are embedded as a one-hot encoding of celltype labels, with a small"
        " amount of random noise added to the embedding"
    ),
)
def celltype_random_embedding_jitter(adata, test=False):
    adata.obsm["X_emb"] = _random_embedding(partition=adata.obs["labels"], jitter=0.01)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Embedding by Celltype",
    method_summary="Cells are embedded as a one-hot encoding of celltype labels",
)
def celltype_random_embedding(adata, test=False):
    adata.obsm["X_emb"] = _random_embedding(partition=adata.obs["labels"], jitter=None)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="No Integration by Batch",
    method_summary="Cells are embedded by computing PCA independently on each batch",
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

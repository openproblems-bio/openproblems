from .....tools.decorators import method
from .....tools.utils import check_version
from ..._common.methods.baseline import _random_embedding

import functools
import numpy as np
import scanpy as sc

_baseline_method = functools.partial(
    method,
    paper_name="Open Problems for Single Cell Analysis",
    paper_reference="openproblems",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)


@_baseline_method(
    method_name="Random Embedding by Celltype",
)
def celltype_random_embedding(adata, test=False):
    adata.obsm["X_emb"] = _random_embedding(partition=adata.obs["labels"], jitter=None)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@_baseline_method(
    method_name="No Integration by Batch",
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

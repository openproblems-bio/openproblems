from ....tools.decorators import method
from ....tools.utils import check_version
from typing import Optional

import numpy as np


@method(
    method_name="Random Features",
    paper_name="Random Features (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def random_features(adata, test=False):
    adata.obsm["X_emb"] = np.random.normal(0, 1, (adata.shape[0], 2))
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="High-dimensional PCA",
    paper_name="High-dimensional PCA (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def high_dim_pca(adata, n_comps: Optional[int] = None, test=False):
    # We wanted to use all features, but output must be dense
    # so this is a close approximation
    import scanpy as sc

    if test:
        n_comps = n_comps or 50
    else:  # pragma: nocover
        n_comps = n_comps or 500

    sc.pp.pca(adata, n_comps=min(min(adata.shape), n_comps))
    adata.obsm["X_emb"] = adata.obsm["X_pca"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

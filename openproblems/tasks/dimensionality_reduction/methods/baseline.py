from ....tools.decorators import method
from ....tools.normalize import log_cp10k
from ....tools.utils import check_version
from .diffusion_map import diffusion_map
from typing import Optional

import functools
import numpy as np

_baseline_method = functools.partial(
    method,
    paper_name="Open Problems for Single Cell Analysis",
    paper_reference="openproblems",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)


@_baseline_method(
    method_name="Random Features",
)
def random_features(adata, test=False):
    adata.obsm["X_emb"] = np.random.normal(0, 1, (adata.shape[0], 2))
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@_baseline_method(
    method_name="True Features",
)
def true_features(adata, test=False):
    adata = log_cp10k(adata)
    adata.obsm["X_emb"] = adata.X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@_baseline_method(
    method_name="Spectral Features",
)
def spectral_features(adata, test=False, n_comps: Optional[int] = None):

    if test:
        n_comps = n_comps or 20
    else:
        n_comps = n_comps or 1000

    n_comps = min(n_comps, min(adata.shape) - 2)

    return diffusion_map(adata, n_comps=n_comps)

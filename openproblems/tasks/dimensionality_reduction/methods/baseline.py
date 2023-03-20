from ....tools.decorators import baseline_method
from ....tools.normalize import log_cp10k
from ....tools.utils import check_version
from .diffusion_map import diffusion_map
from typing import Optional

import numpy as np


@baseline_method(
    method_name="Random Features",
    method_summary=(
        "Randomly generated two-dimensional coordinates from a normal distribution."
    ),
)
def random_features(adata, test=False):
    adata.obsm["X_emb"] = np.random.normal(0, 1, (adata.shape[0], 2))
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="True Features",
    method_summary="Use of the original feature inputs as the 'embedding'.",
)
def true_features(adata, test=False):
    adata.obsm["X_emb"] = adata.X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="True Features (logCP10k)",
    method_summary="Use of the original feature inputs as the 'embedding'.",
)
def true_features_log_cp10k(adata, test=False):
    adata = log_cp10k(adata)
    adata.obsm["X_emb"] = adata.X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Spectral Features",
    method_summary="Use 1000-dimensional diffusions maps as an embedding",
)
def spectral_features(adata, test=False, n_comps: Optional[int] = None):

    if test:
        n_comps = n_comps or 20
    else:
        n_comps = n_comps or 1000

    n_comps = min(n_comps, min(adata.shape) - 2)

    return diffusion_map(adata, n_comps=n_comps)

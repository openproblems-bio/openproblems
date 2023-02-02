from ....tools.decorators import baseline_method
from ....tools.normalize import log_cp10k
from ....tools.normalize import log_cp10k_hvg
from ....tools.utils import check_version

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
    method_name="True Features (logCP10k, 1kHVG)",
    method_summary="Use of the original feature inputs as the 'embedding'.",
)
def true_features_log_cp10k_hvg(adata, test=False):
    adata = log_cp10k_hvg(adata)
    adata.obsm["X_emb"] = adata[:, adata.var["highly_variable"]].copy().X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

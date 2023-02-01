from ....tools.decorators import baseline_method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

import numpy as np


@baseline_method(method_name="Random Features", method_summary="TODO")
def random_features(adata, test=False):
    adata.obsm["X_emb"] = np.random.normal(0, 1, (adata.shape[0], 2))
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(method_name="True Features", method_summary="TODO")
def true_features(adata, test=False):
    adata.obsm["X_emb"] = adata.X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(method_name="True Features (logCPM)", method_summary="TODO")
def true_features_log_cpm(adata, test=False):
    adata = log_cpm(adata)
    adata.obsm["X_emb"] = adata.X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(method_name="True Features (logCPM, 1kHVG)", method_summary="TODO")
def true_features_log_cpm_hvg(adata, test=False):
    adata = log_cpm_hvg(adata)
    adata.obsm["X_emb"] = adata[:, adata.var["highly_variable"]].copy().X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

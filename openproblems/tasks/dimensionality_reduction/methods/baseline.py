from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_cpm_hvg
from ....tools.utils import check_version

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
    adata.obsm["X_emb"] = adata.X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@_baseline_method(
    method_name="True Features (logCPM)",
)
def true_features_log_cpm(adata, test=False):
    adata = log_cpm(adata)
    adata.obsm["X_emb"] = adata.X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@_baseline_method(
    method_name="True Features (logCPM, 1kHVG)",
)
def true_features_log_cpm_hvg(adata, test=False):
    adata = log_cpm_hvg(adata)
    adata = adata[:, adata.var["highly_variable"]].copy()
    adata.obsm["X_emb"] = adata.X
    if test:
        adata.obsm["X_emb"] = adata.obsm["X_emb"][:, :100]

    adata.obsm["X_emb"] = adata.obsm["X_emb"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

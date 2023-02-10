from .....tools.decorators import method
from .....tools.utils import check_version
from ...batch_integration_graph.methods.baseline import _randomize_features
from ...batch_integration_graph.methods.baseline import _randomize_graph
from ...batch_integration_graph.methods.baseline import _set_uns

import functools

_baseline_method = functools.partial(
    method,
    paper_name="Open Problems for Single Cell Analysis",
    paper_reference="openproblems",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)


@_baseline_method(
    method_name="No Integration",
)
def no_integration(adata, test=False):
    adata.obsp["connectivities"] = adata.obsp["uni_connectivities"]
    adata.obsp["distances"] = adata.obsp["uni_distances"]
    _set_uns(adata)
    adata.obsm["X_emb"] = adata.obsm["X_uni_pca"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@_baseline_method(
    method_name="Random Integration",
)
def random_integration(adata, test=False):
    adata.X = _randomize_features(adata.X)
    adata.obsm["X_emb"] = _randomize_features(adata.obsm["X_uni_pca"])
    adata = _randomize_graph(adata)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@_baseline_method(
    method_name="Random Integration by Celltype",
    paper_name="Random Integration by Celltype (baseline)",
    paper_reference="openproblems",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    is_baseline=True,
)
def celltype_random_integration(adata, test=False):
    adata.obsm["X_emb"] = _randomize_features(
        adata.obsm["X_uni_pca"], partition=adata.obs["labels"]
    )
    adata.X = _randomize_features(adata.X, partition=adata.obs["labels"])
    adata = _randomize_graph(
        adata,
        partition=adata.obs["labels"].to_numpy(),
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@_baseline_method(
    method_name="Random Integration by Batch",
)
def batch_random_integration(adata, test=False):
    adata.obsm["X_emb"] = _randomize_features(
        adata.obsm["X_uni_pca"], partition=adata.obs["batch"]
    )
    adata.X = _randomize_features(adata.X, partition=adata.obs["batch"])
    adata = _randomize_graph(
        adata,
        partition=adata.obs["batch"].to_numpy(),
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

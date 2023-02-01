from .....tools.decorators import baseline_method
from .....tools.utils import check_version
from ...batch_integration_graph.methods.baseline import _randomize_features


@baseline_method(
    method_name="No Integration",
    method_summary="TODO",
)
def no_integration(adata, test=False):
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Integration",
    method_summary="TODO",
)
def random_integration(adata, test=False):
    adata.X = _randomize_features(adata.X)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Integration by Celltype",
    method_summary="TODO",
)
def celltype_random_integration(adata, test=False):
    adata.X = _randomize_features(adata.X, partition=adata.obs["labels"])
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Integration by Batch",
    method_summary="TODO",
)
def batch_random_integration(adata, test=False):
    adata.X = _randomize_features(adata.X, partition=adata.obs["batch"])
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

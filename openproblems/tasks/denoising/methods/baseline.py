from ....tools.decorators import baseline_method
from ....tools.utils import check_version


@baseline_method(
    method_name="No denoising",
    method_summary="Denoised outputs are defined from the unmodified input data.",
)
def no_denoising(adata, test=False):
    """Do nothing."""
    adata.obsm["denoised"] = adata.obsm["train"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Perfect denoising",
    method_summary="Denoised outputs are defined from the target data.",
)
def perfect_denoising(adata, test=False):
    """Cheat."""
    adata.obsm["denoised"] = adata.obsm["test"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

from ....tools.decorators import method
from ....tools.utils import check_version


@method(
    method_name="No denoising",
    paper_name="Molecular Cross-Validation for Single-Cell RNA-seq",
    paper_reference="https://doi.org/10.1101/786269",
    paper_year=2019,
    code_url="https://github.com/czbiohub/molecular-cross-validation",
    is_baseline=True,
)
def no_denoising(adata, test=False):
    """Do nothing."""
    adata.obsm["denoised"] = adata.obsm["train"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@method(
    method_name="Perfect denoising",
    paper_name="Molecular Cross-Validation for Single-Cell RNA-seq",
    paper_reference="https://doi.org/10.1101/786269",
    paper_year=2019,
    code_url="https://github.com/czbiohub/molecular-cross-validation",
    is_baseline=True,
)
def perfect_denoising(adata, test=False):
    """Cheat."""
    adata.obsm["denoised"] = adata.obsm["test"].toarray()
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata

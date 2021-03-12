# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version

# We define one method for each scenario (scaling/hvg)
# HVG Selection with `_hvg.hvg_batch` to not subset in testcases,
# wrapper around `scIB.pp.hvg_batch()`
# Scaling with `scIB.pp.scale_batch`


@method(
    method_name="Example (full, unscaled)",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="",
    code_version=check_version("bbknn"),
    image="openproblems-python-batch-integration",  # only if required
)
def _example_full_unscaled(adata):
    from scIB.integration import runBBKNN

    adata = runBBKNN(adata, "batch")
    # Complete the result in-place
    return adata


@method(
    method_name="Example (hvg,unscaled)",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="",
    code_version=check_version("bbknn"),
    image="openproblems-python-batch-integration",  # only if required
)
def _example_hvg_unscaled(adata):
    from _hvg import hvg_batch
    from scIB.integration import runBBKNN

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runBBKNN(adata, "batch")
    # Complete the result in-place
    return adata


@method(
    method_name="Example (hvg,scaled)",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="",
    code_version=check_version("bbknn"),
    image="openproblems-python-batch-integration",  # only if required
)
def _example_hvg_scaled(adata):
    from _hvg import hvg_batch
    from scIB.integration import runBBKNN
    from scIB.preprocessing import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runBBKNN(adata, "batch")
    # Complete the result in-place
    return adata


@method(
    method_name="Example (full/scaled)",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="",
    code_version=check_version("bbknn"),
    image="openproblems-python-batch-integration",  # only if required
)
def _example_full_scaled(adata):
    from scIB.integration import runBBKNN
    from scIB.preprocessing import scale_batch

    adata = scale_batch(adata, "batch")
    adata = runBBKNN(adata, "batch")
    # Complete the result in-place
    return adata

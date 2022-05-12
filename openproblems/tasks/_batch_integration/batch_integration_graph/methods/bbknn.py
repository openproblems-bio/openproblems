# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version

def _run_bbknn(adata, batch):
    import bbknn
    from scanpy.preprocessing import pca

    pca(adata, svd_solver="arpack")
    if adata.n_obs < 1e5:
        return bbknn.bbknn(adata, batch_key=batch, copy=True)
    if adata.n_obs >= 1e5:
        return bbknn.bbknn(
            adata, batch_key=batch, neighbors_within_batch=25, copy=True
        )



@method(
    method_name="BBKNN",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="https://github.com/Teichlab/bbknn",
    code_version=check_version("bbknn"),
    image="openproblems-python-batch-integration",  # only if required
)
def bbknn_full_unscaled(adata, test=False):

    adata = _run_bbknn(adata, "batch")
    # Complete the result in-place
    return adata


@method(
    method_name="BBKNN (hvg,unscaled)",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="https://github.com/Teichlab/bbknn",
    code_version=check_version("bbknn"),
    image="openproblems-python-batch-integration",  # only if required
)
def bbknn_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = _run_bbknn(adata, "batch")
    # Complete the result in-place
    return adata


@method(
    method_name="BBKNN (hvg,scaled)",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="https://github.com/Teichlab/bbknn",
    code_version=check_version("bbknn"),
    image="openproblems-python-batch-integration",  # only if required
)
def bbknn_hvg_scaled(adata, test=False):
    from ._utils import hvg_batch
    from ._utils import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = _run_bbknn(adata, "batch")
    # Complete the result in-place
    return adata


@method(
    method_name="BBKNN (full/scaled)",
    paper_name="BBKNN: fast batch alignment of single cell transcriptomes",
    paper_url="https://academic.oup.com/bioinformatics/article/36/3/964/5545955",
    paper_year=2020,
    code_url="https://github.com/Teichlab/bbknn",
    code_version=check_version("bbknn"),
    image="openproblems-python-batch-integration",  # only if required
)
def bbknn_full_scaled(adata, test=False):
    from ._utils import scale_batch

    adata = scale_batch(adata, "batch")
    adata = _run_bbknn(adata, "batch")
    # Complete the result in-place
    return adata

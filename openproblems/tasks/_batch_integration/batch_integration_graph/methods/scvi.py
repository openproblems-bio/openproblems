from .....tools.decorators import method
from .....tools.utils import check_version
from typing import Optional

import functools

_scvi_method = functools.partial(
    method,
    paper_name="Deep generative modeling for single-cell transcriptomics",
    paper_url="https://www.nature.com/articles/s41592-018-0229-2",
    paper_year=2018,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-python-batch-integration",  # only if required
)


def _scvi(adata, test: bool = False, max_epochs: Optional[int] = None):
    from openproblems.tools.normalize import log_scran_pooling
    from scanpy.preprocessing import neighbors
    from scib.integration import scvi

    if test:
        max_epochs = max_epochs or 2

    adata = log_scran_pooling(adata)

    adata.obs.rename(
        columns={"labels": "lab"}, inplace=True
    )  # ugly fix for scvi conversion error
    adata = scvi(adata, "batch", max_epochs=max_epochs)
    neighbors(adata, use_rep="X_emb")
    adata.obs.rename(
        columns={"lab": "labels"}, inplace=True
    )  # ugly fix for scvi conversion error
    # Complete the result in-place
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata


@_scvi_method(method_name="scVI (full/unscaled)")
def scvi_full_unscaled(adata, test: bool = False, max_epochs: Optional[int] = None):
    return _scvi(adata, test=test, max_epochs=max_epochs)


@_scvi_method(method_name="scVI (hvg/unscaled)")
def scvi_hvg_unscaled(adata, test: bool = False, max_epochs: Optional[int] = None):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _scvi(adata, test=test, max_epochs=max_epochs)

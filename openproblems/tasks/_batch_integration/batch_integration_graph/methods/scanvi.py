from .....tools.decorators import method
from .....tools.utils import check_version
from typing import Optional

import functools

_scanvi_method = functools.partial(
    method,
    paper_name="Probabilistic harmonization and annotation of single‐cell "
    "transcriptomics data with deep generative models",
    paper_url="https://doi.org/10.15252/msb.20209620",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-r-pytorch",
)


def _scanvi(adata, test: bool = False, max_epochs: Optional[int] = None):
    from scanpy.preprocessing import neighbors
    from scib.integration import scanvi

    if test:
        max_epochs = max_epochs or 2

    adata.obs.rename(
        columns={"labels": "lab"}, inplace=True
    )  # ugly fix for scvi conversion error
    adata = scanvi(adata, "batch", "lab", max_epochs=max_epochs)
    neighbors(adata, use_rep="X_emb")
    adata.obs.rename(
        columns={"lab": "labels"}, inplace=True
    )  # ugly fix for scvi conversion error
    adata.uns["method_code_version"] = check_version("scvi-tools")
    return adata


@_scanvi_method(method_name="scANVI (full/unscaled)")
def scanvi_full_unscaled(adata, test: bool = False, max_epochs: Optional[int] = None):
    return _scanvi(adata, test=test, max_epochs=max_epochs)


@_scanvi_method(method_name="scANVI (hvg/unscaled)")
def scanvi_hvg_unscaled(adata, test: bool = False, max_epochs: Optional[int] = None):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _scanvi(adata, test=test, max_epochs=max_epochs)

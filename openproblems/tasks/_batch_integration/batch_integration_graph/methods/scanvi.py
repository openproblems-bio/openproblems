from .....tools.decorators import method
from .....tools.utils import check_version

import functools

_scanvi_method = functools.partial(
    method,
    paper_name="Probabilistic harmonization and annotation of single‐cell "
    "transcriptomics data with deep generative models",
    paper_url="https://www.embopress.org/doi/full/10.15252/msb.20209620",
    paper_year=2021,
    code_url="https://github.com/YosefLab/scvi-tools",
    image="openproblems-python-batch-integration",
)


def _scanvi(adata):
    from scanpy.preprocessing import neighbors
    from scib.integration import runScanvi

    adata.obs.rename(
        columns={"labels": "lab"}, inplace=True
    )  # ugly fix for scvi conversion error
    adata = runScanvi(adata, "batch", "lab")
    neighbors(adata, use_rep="X_emb")
    adata.obs.rename(
        columns={"lab": "labels"}, inplace=True
    )  # ugly fix for scvi conversion error
    adata.uns["method_code_version"] = check_version("scvi")
    return adata


@_scanvi_method(method_name="scANVI (full/unscaled)")
def scanvi_full_unscaled(adata, test=False):
    return _scanvi(adata)


@_scanvi_method(method_name="scANVI (hvg/unscaled)")
def scanvi_hvg_unscaled(adata, test=False):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _scanvi(adata)

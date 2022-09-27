from .....tools.conversion import r_function
from .....tools.decorators import method
from .....tools.utils import check_r_version
from typing import Optional

import functools

_liger_method = functools.partial(
    method,
    paper_name="Single-Cell Multi-omic Integration Compares and "
    "Contrasts Features of Brain Cell Identity",
    paper_url="https://doi.org/10.1016/j.cell.2019.05.006",
    paper_year=2019,
    code_url="https://github.com/welch-lab/liger",
    image="openproblems-r-extras",
)

_r_liger = r_function("liger.R", args="sce, batch, k, nrep, thresh")


def _liger(
    adata,
    batch: str,
    test: bool = False,
    k: Optional[int] = None,
    nrep: Optional[int] = None,
    thresh: Optional[float] = None,
):
    from openproblems.tools.normalize import log_scran_pooling
    from scanpy.pp import neighbors

    if test:
        k = k or 5
        nrep = nrep or 1
        thresh = thresh or 5e-3
    else:  # pragma: nocover
        k = k or 20
        nrep = nrep or 3
        thresh = thresh or 5e-5

    adata.obsm["X_emb"] = _r_liger(adata, batch, k, nrep, thresh)
    neighbors(adata, use_rep="X_emb")

    adata.uns["method_code_version"] = check_r_version("liger")
    return adata


@_liger_method(method_name="Liger (full/unscaled)")
def liger_full_unscaled(
    adata,
    test: bool = False,
    k: Optional[int] = None,
    nrep: Optional[int] = None,
    thresh: Optional[float] = None,
):
    return _liger(adata, "batch", test=test, k=k, nrep=nrep, thresh=thresh)


@_liger_method(method_name="Liger (hvg/unscaled)")
def liger_hvg_unscaled(
    adata,
    test: bool = False,
    k: Optional[int] = None,
    nrep: Optional[int] = None,
    thresh: Optional[float] = None,
):
    from ._utils import hvg_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    return _liger(adata, "batch", test=test, k=k, nrep=nrep, thresh=thresh)

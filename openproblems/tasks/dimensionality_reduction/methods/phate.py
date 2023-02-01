from ....tools.decorators import method
from ....tools.normalize import log_cpm
from ....tools.normalize import log_cpm_hvg
from ....tools.normalize import sqrt_cpm
from ....tools.utils import check_version
from typing import Optional

import functools

_phate_method = functools.partial(
    method,
    method_summary="TODO",
    paper_name=(
        "Visualizing Structure and Transitions in High-Dimensional Biological Data"
    ),
    paper_reference="moon2019visualizing",
    paper_year=2019,
    code_url="https://github.com/KrishnaswamyLab/PHATE/",
    image="openproblems-python-extras",
)


def _phate(
    adata, test: bool = False, genes=None, n_pca: Optional[int] = None, gamma: float = 1
):
    from phate import PHATE

    if test:
        n_pca = n_pca or 10
    else:  # pragma: no cover
        n_pca = n_pca or 100

    if genes is not None:
        X = adata[:, genes].copy().X
    else:
        X = adata.X

    phate_op = PHATE(n_pca=n_pca, verbose=False, n_jobs=-1, gamma=gamma)
    adata.obsm["X_emb"] = phate_op.fit_transform(X)
    adata.uns["method_code_version"] = check_version("phate")
    return adata


@_phate_method(method_name="PHATE (default)")
def phate_default(adata, test: bool = False, n_pca: Optional[int] = None):
    adata = sqrt_cpm(adata)
    adata = _phate(adata, test=test, n_pca=n_pca)
    # revert to expected adata.X
    adata = log_cpm(adata)
    return adata


@_phate_method(method_name="PHATE (gamma=0)")
def phate_sqrt(adata, test: bool = False, n_pca: Optional[int] = None):
    adata = sqrt_cpm(adata)
    adata = _phate(adata, test=test, n_pca=n_pca, gamma=0)
    # revert to expected adata.X
    adata = log_cpm(adata)
    return adata


@_phate_method(method_name="PHATE (logCPM)")
def phate_logCPM_1kHVG(adata, test: bool = False, n_pca: Optional[int] = None):
    adata = log_cpm(adata)
    return _phate(adata, test=test, n_pca=n_pca)


@_phate_method(method_name="PHATE (logCPM, 1kHVG)")
def phate_logCPM(adata, test: bool = False, n_pca: Optional[int] = None):
    adata = log_cpm_hvg(adata)
    return _phate(adata, test=test, genes=adata.var["highly_variable"], n_pca=n_pca)

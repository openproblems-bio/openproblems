from ....tools.decorators import method
from ....tools.normalize import log_cp10k
from ....tools.normalize import log_cp10k_hvg
from ....tools.normalize import sqrt_cp10k
from ....tools.utils import check_version
from typing import Optional

import functools

_phate_method = functools.partial(
    method,
    method_summary=(
        "PHATE or “Potential of Heat - diffusion for Affinity - based Transition"
        " Embedding” uses the potential of heat diffusion to preserve trajectories in a"
        " dataset via a diffusion process. It is an affinity - based method that"
        " creates an embedding by finding the dominant eigenvalues of a Markov"
        " transition matrix. We evaluate several variants including using the"
        " recommended square - root transformed CPM matrix as input, this input with"
        " the gamma parameter set to zero and the normal logCPM transformed matrix with"
        " and without HVG selection."
    ),
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
    adata = sqrt_cp10k(adata)
    adata = _phate(adata, test=test, n_pca=n_pca)
    # revert to expected adata.X
    adata = log_cp10k(adata)
    return adata


@_phate_method(method_name="PHATE (gamma=0)")
def phate_sqrt(adata, test: bool = False, n_pca: Optional[int] = None):
    adata = sqrt_cp10k(adata)
    adata = _phate(adata, test=test, n_pca=n_pca, gamma=0)
    # revert to expected adata.X
    adata = log_cp10k(adata)
    return adata


@_phate_method(method_name="PHATE (logCP10k)")
def phate_logCP10k_1kHVG(adata, test: bool = False, n_pca: Optional[int] = None):
    adata = log_cp10k(adata)
    return _phate(adata, test=test, n_pca=n_pca)


@_phate_method(method_name="PHATE (logCP10k, 1kHVG)")
def phate_logCP10k(adata, test: bool = False, n_pca: Optional[int] = None):
    adata = log_cp10k_hvg(adata)
    return _phate(adata, test=test, genes=adata.var["highly_variable"], n_pca=n_pca)

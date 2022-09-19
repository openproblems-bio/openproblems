from ....tools.decorators import method
from ....tools.utils import check_version

import functools
import numpy as np
import scprep

_magic_method = functools.partial(
    method,
    paper_name="Recovering Gene Interactions from Single-Cell Data "
    "Using Data Diffusion",
    paper_url="https://doi.org/10.1016/j.cell.2018.05.061",
    paper_year=2018,
    code_url="https://github.com/KrishnaswamyLab/MAGIC",
    image="openproblems-python-extras",
)


def _magic(adata, solver, normtype, decay, t):
    from magic import MAGIC

    X, libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train"], rescale=1, return_library_size=True
    )

    if normtype == "sqrt":
        X = scprep.transform.sqrt(X)
    elif normtype == "log":
        X = scprep.transform.log(X, base="e")
    Y = MAGIC(solver=solver, decay=decay, t=t, verbose=False).fit_transform(
        X, genes="all_genes"
    )
    if normtype == "sqrt":
        Y = scprep.utils.matrix_transform(Y, np.square)
    elif normtype == "log":
        Y = scprep.utils.matrix_transform(Y, np.expm1)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y

    adata.uns["method_code_version"] = check_version("magic-impute")
    return adata


@_magic_method(
    method_name="MAGIC",
)
def magic(adata, test=False, normtype="sqrt", decay=1, t=3):
    return _magic(adata, solver="exact")


@method(
    method_name="KNN Smoothing",
    paper_name="KNN Smoothing (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
)
def knn_naive(adata, test=False, normtype="log", decay=0, t=1):
    return _magic(adata, solver="exact")


@_magic_method(
    method_name="MAGIC (approximate)",
)
def magic_approx(adata, test=False):
    return _magic(adata, solver="approximate")

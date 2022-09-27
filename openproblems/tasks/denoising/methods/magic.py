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


def _magic(adata, solver, normtype="sqrt", **kwargs):
    from magic import MAGIC

    if normtype == "sqrt":
        norm_fn = np.sqrt
        denorm_fn = np.square
    elif normtype == "log":
        norm_fn = np.log1p
        denorm_fn = np.expm1
    else:
        raise NotImplementedError

    X, libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train"], rescale=1, return_library_size=True
    )

    X = scprep.utils.matrix_transform(X, norm_fn)
    Y = MAGIC(solver=solver, **kwargs, verbose=False).fit_transform(
        X, genes="all_genes"
    )

    Y = scprep.utils.matrix_transform(Y, denorm_fn)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y

    adata.uns["method_code_version"] = check_version("magic-impute")
    return adata


@_magic_method(
    method_name="MAGIC",
)
def magic(adata, test=False):
    return _magic(adata, solver="exact", normtype="sqrt")


@_magic_method(
    method_name="MAGIC (approximate)",
)
def magic_approx(adata, test=False):
    return _magic(adata, solver="approximate", normtype="sqrt")


@method(
    method_name="KNN Smoothing",
    paper_name="KNN Smoothing (baseline)",
    paper_url="https://openproblems.bio",
    paper_year=2022,
    code_url="https://github.com/openproblems-bio/openproblems",
    image="openproblems-python-extras",
)
def knn_naive(adata, test=False):
    return _magic(adata, solver="exact", normtype="log", decay=None, t=1)

from ....tools.decorators import baseline_method
from ....tools.decorators import method
from ....tools.utils import check_version

import functools
import numpy as np
import scprep

_magic_method = functools.partial(
    method,
    method_summary="TODO",
    paper_name=(
        "Recovering Gene Interactions from Single-Cell Data Using Data Diffusion"
    ),
    paper_reference="van2018recovering",
    paper_year=2018,
    code_url="https://github.com/KrishnaswamyLab/MAGIC",
    image="openproblems-python-extras",
)


def _magic(adata, solver, normtype="sqrt", reverse_norm_order=False, **kwargs):
    from magic import MAGIC

    if normtype == "sqrt":
        norm_fn = np.sqrt
        denorm_fn = np.square
    elif normtype == "log":
        norm_fn = np.log1p
        denorm_fn = np.expm1
    else:
        raise NotImplementedError

    X = adata.obsm["train"]
    if reverse_norm_order:
        # inexplicably, this sometimes performs better
        X = scprep.utils.matrix_transform(X, norm_fn)
        X, libsize = scprep.normalize.library_size_normalize(
            X, rescale=1, return_library_size=True
        )
    else:
        X, libsize = scprep.normalize.library_size_normalize(
            X, rescale=1, return_library_size=True
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
    method_name="MAGIC (reversed normalization)",
)
def magic_reverse_norm(adata, test=False):
    return _magic(adata, solver="exact", normtype="sqrt", reverse_norm_order=True)


@_magic_method(
    method_name="MAGIC (approximate)",
)
def magic_approx(adata, test=False):
    return _magic(adata, solver="approximate", normtype="sqrt")


@_magic_method(
    method_name="MAGIC (approximate, reversed normalization)",
)
def magic_approx_reverse_norm(adata, test=False):
    return _magic(adata, solver="approximate", normtype="sqrt", reverse_norm_order=True)


@baseline_method(
    method_name="KNN smoothing",
    method_summary="TODO",
    is_baseline=False,
    image="openproblems-python-extras",
)
def knn_naive(adata, test=False):
    return _magic(adata, solver="exact", normtype="log", decay=None, t=1)

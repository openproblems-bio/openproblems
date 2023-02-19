from ....tools.conversion import r_function
from ....tools.decorators import method

import logging

_alra = r_function("alra.R")

log = logging.getLogger("openproblems")


method_name = ("ALRA (sqrt norm, reversed normalization)",)
_alra_method = functools.partial(
    method,
    paper_name="Zero-preserving imputation of scRNA-seq data using "
    "low-rank approximation",
    paper_reference="linderman2018zero",
    paper_year=2018,
    code_url="https://github.com/KlugerLab/ALRA",
    image="openproblems-r-extras",
)


def _alra(adata, normtype="log", reverse_norm_order=False, test=False):
    import numpy as np
    import rpy2.rinterface_lib.embedded
    import scprep

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
        adata.obsm["train"] = scprep.utils.matrix_transform(
            adata.obsm["train"], norm_fn
        )
        adata.obsm["train"], libsize = scprep.normalize.library_size_normalize(
            adata.obsm["train"], rescale=1, return_library_size=True
        )
    else:
        adata.obsm["train"], libsize = scprep.normalize.library_size_normalize(
            adata.obsm["train"], rescale=1, return_library_size=True
        )
        adata.obsm["train"] = scprep.utils.matrix_transform(X, norm_fn)

    adata.obsm["train"] = adata.obsm["train"].tocsr()
    # run alra
    # _alra takes sparse array, returns dense array
    Y = None
    attempts = 0
    while Y is None:
        try:
            Y = _alra(adata)
        except rpy2.rinterface_lib.embedded.RRuntimeError:  # pragma: no cover
            if attempts < 10:
                attempts += 1
                log.warning(f"alra.R failed (attempt {attempts})")
            else:
                raise

    # transform back into original space
    # functions are reversed!
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y

    adata.uns["method_code_version"] = "1.0.0"
    return adata


@_alra_method(
    method_name="ALRA (sqrt norm, reversed normalization)",
)
def alra_sqrt_reversenorm(adata, test=False):
    return _alra(adata, normtype="log", reverse_norm_order=True, test=False)


@_alra_method(
    method_name="ALRA (log norm, reversed normalization)",
)
def alra_log_reversenorm(adata, test=False):
    return _alra(adata, normtype="log", reverse_norm_order=True, test=False)


@_alra_method(
    method_name="ALRA (sqrt norm)",
)
def alra_sqrt(adata, test=False):
    return _alra(adata, normtype="log", reverse_norm_order=False, test=False)


@_alra_method(
    method_name="ALRA (log norm)",
)
def alra_log(adata, test=False):
    return _alra(adata, normtype="log", reverse_norm_order=False, test=False)

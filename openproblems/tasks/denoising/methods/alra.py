from ....tools.conversion import r_function
from ....tools.decorators import method

_alra = r_function("alra.R")


@method(
    method_name="ALRA",
    paper_name="Zero-preserving imputation of scRNA-seq...",
    paper_url="https://www.biorxiv.org/content/10.1101/397588v1",
    paper_year=2018,
    code_url="https://github.com/KlugerLab/ALRA",
    code_version="v1.0.0",
    image="openproblems-r-extras",
)
def alra(adata, test=False):
    import numpy as np
    import rpy2.rinterface_lib.embedded
    import scprep

    # libsize and sqrt norm
    adata.obsm["train_norm"] = scprep.utils.matrix_transform(
        adata.obsm["train"], np.sqrt
    )
    adata.obsm["train_norm"], libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train_norm"], rescale=1, return_library_size=True
    )
    # run alra
    # _alra takes sparse array, returns dense array
    Y = None
    attempts = 0
    while Y is None:
        try:
            Y = _alra(adata.obsm["train_norm"])
        except rpy2.rinterface_lib.embedded.RRuntimeError as e:
            if "non-comfortable arguments" in str(e) and attempts < 5:
                attempts += 1
            else:
                raise

    # transform back into original space
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y
    return adata

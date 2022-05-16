from ....tools.conversion import r_function
from ....tools.decorators import method

import logging

_alra = r_function("alra.R")

log = logging.getLogger("openproblems")


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
    data_norm = scprep.utils.matrix_transform(adata.obsm["train"], np.sqrt)
    data_norm, libsize = scprep.normalize.library_size_normalize(
        data_norm, rescale=1, return_library_size=True
    )
    data_norm = data_norm.tocsc()
    # run alra
    # _alra takes sparse array, returns dense array
    Y = None
    attempts = 0
    while Y is None:
        try:
            Y = _alra(data_norm)
        except rpy2.rinterface_lib.embedded.RRuntimeError as e:
            if attempts < 5:
                log.warning("alra.R failed")
                log.warning(str(e))
                attempts += 1
            else:
                raise

    # transform back into original space
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y
    return adata

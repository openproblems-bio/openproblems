from ....tools.conversion import r_function
from ....tools.decorators import method

# should be imported by docker
import numpy as np
import scprep

# from ....tools.utils import check_version


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
    # libsize and sqrt norm
    # overwrite obsm['train'] with normalized version
    adata.obsm["train"] = scprep.utils.matrix_transform(adata.obsm["train"], np.sqrt)
    adata.obsm["train"], libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train"], rescale=1, return_library_size=True
    )
    # run alra
    # _alra takes adata returns adata, edits "train"
    Y = _alra(adata)
    Y = Y.obsm["train"]
    # transform back into original space

    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y
    return adata

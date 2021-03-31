from ....tools.conversion import r_function
from ....tools.decorators import method
from ....tools.utils import check_version

# should be imported by docker
import numpy as np
import scprep

_alra = r_function("alra.R")


@method(
    method_name="ALRA",
    paper_name="Zero-preserving imputation of scRNA-seq...",
    paper_url="https://www.biorxiv.org/content/10.1101/397588v1",
    paper_year=2018,
    code_url="https://github.com/KlugerLab/ALRA",
    code_version=check_version("scprep"),
    image="openproblems-r-extras",
)
def alra(adata):
    # libsize and sqrt norm
    # overwrite obsm['train'] with normalized version
    adata.obsm["train"], libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train"], rescale=1, return_library_size=True
    )

    Y = _alra(adata)["train"]
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y
    return adata

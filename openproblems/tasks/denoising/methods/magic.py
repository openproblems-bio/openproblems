from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np
import scprep


@method(
    method_name="MAGIC",
    paper_name="Recovering Gene Interactions from Single-Cell Data "
    "Using Data Diffusion",
    paper_url="https://www.cell.com/cell/abstract/S0092-8674(18)30724-4",
    paper_year=2018,
    code_url="https://github.com/KrishnaswamyLab/MAGIC",
    code_version=check_version("magic-impute"),
    image="openproblems-python-extras",
)
def magic(adata):
    from magic import MAGIC

    X, libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train"], rescale=1, return_library_size=True
    )
    X = scprep.transform.sqrt(X)
    Y = MAGIC().fit_transform(X, genes="all_genes")
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y
    return adata

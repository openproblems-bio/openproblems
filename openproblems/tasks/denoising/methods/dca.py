from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np
import scprep


def _dca(adata):
    from dca.api import dca

    X, libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train"], rescale=1, return_library_size=True
    )
    X = scprep.transform.sqrt(X)
    Y = dca(adata, threads=1)
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y
    adata.uns["method_code_version"] = check_version("dca")
    return adata


@method(
    method_name="DCA",
    paper_name="Single-cell RNA-seq denoising using...",
    paper_url="https://www.nature.com/articles/s41467-018-07931-2",
    paper_year=2019,
    code_url="https://github.com/theislab/dca",
    image="openproblems-python-tf2.4",
)
def dca(adata, test=False):
    return _dca(adata)

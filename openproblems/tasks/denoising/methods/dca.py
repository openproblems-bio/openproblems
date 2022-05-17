from ....tools.decorators import method
from ....tools.utils import check_version

import numpy as np
import scprep
import dca.api


def _dca(adata):
    from magic import MAGIC

    X, libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train"], rescale=1, return_library_size=True
    )
    X = scprep.transform.sqrt(X)
    Y = dca(adata_ae, threads=1)
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y
    return adata


@method(
    method_name="DCA",
    paper_name="Single-cell RNA-seq denoising using a deep count autoencoder"
    paper_url="https://www.nature.com/articles/s41467-018-07931-2",
    paper_year=2019,
    code_url="https://github.com/theislab/dca",
    code_version=check_version("dca.api"),
    image="openproblems-python-extras",
)
def magic(adata, test=False):
    return _magic(adata, solver="exact")


@method(
    method_name="MAGIC (approximate)",
    paper_name="Recovering Gene Interactions from Single-Cell Data "
    "Using Data Diffusion",
    paper_url="https://www.cell.com/cell/abstract/S0092-8674(18)30724-4",
    paper_year=2018,
    code_url="https://github.com/KrishnaswamyLab/MAGIC",
    code_version=check_version("magic-impute"),
    image="openproblems-python-extras",
)
def dca_run(adata, test=False):
    return _dca(adata)

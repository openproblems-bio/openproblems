from ....tools.decorators import method
from ....tools.utils import check_version

import functools
import numpy as np
import scprep


def _magic(adata, solver):
    from magic import MAGIC

    X, libsize = scprep.normalize.library_size_normalize(
        adata.obsm["train"], rescale=1, return_library_size=True
    )
    X = scprep.transform.sqrt(X)
    #the decay and t parameter give knn smoothing behavior
    Y = MAGIC(solver=solver, decay=0, 
              t=1, verbose=False).fit_transform(X, genes="all_genes")
    Y = scprep.utils.matrix_transform(Y, np.square)
    Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)
    adata.obsm["denoised"] = Y

    adata.uns["method_code_version"] = check_version("magic-impute")
    return adata


def knn_smoothing(adata, test=False):
    return _magic(adata, solver="exact")


@method(
    method_name="KNN smoothing",
    paper_name="K-nearest neighbor smoothing for high-throughput "
    "single-cell RNA-Seq data",
    paper_url="https://www.biorxiv.org/content/10.1101/217737v3",
    paper_year=2018,
    #magic codebase is used instead of knn-smoothing
    code_url="https://github.com/yanailab/knn-smoothing",
    image="openproblems-python-extras",
)

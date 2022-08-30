from ....tools.decorators import method
from ....tools.utils import check_version

import functools
import numpy as np
import scprep


# the decay and t parameter give knn smoothing behavior
def _knn_smoothing(adata, solver):
    import knn-smooth
    
    smoothX = knn-smooth.knn_smoothing(adata.X)
    adata.X = smoothX
    return adata


def knn_smoothing(adata, test=False):
    return _knn_smoothing(adata)

# magic codebase is used instead of knn-smoothing


@method(
    method_name="KNN smoothing",
    paper_name="K-nearest neighbor smoothing for high-throughput "
    "single-cell RNA-Seq data",
    paper_url="https://www.biorxiv.org/content/10.1101/217737v3",
    paper_year=2018,
    code_url="https://github.com/yanailab/knn-smoothing",
    image="openproblems-python-extras")

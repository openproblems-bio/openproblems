from ....tools.decorators import method
from ....tools.utils import check_version

import functools
import numpy as np
import scprep

@method(
    method_name="KNN smoothing",
    paper_name="K-nearest neighbor smoothing for high-throughput "
    "single-cell RNA-Seq data",
    paper_url="https://www.biorxiv.org/content/10.1101/217737v3",
    paper_year=2018,
    code_url="https://github.com/yanailab/knn-smoothing",
    image="openproblems-python-extras")


def knn_smoothing(adata, test=False):
    import knn_smooth
    adata.X = knn_smooth.knn_smoothing(adata.X)
    return adata

from ....tools.decorators import method
from ....tools.utils import check_version
import numpy as np


@method(
    method_name="KNN smoothing",
    paper_name="K-nearest neighbor smoothing for high-throughput "
    "single-cell RNA-Seq data",
    paper_url="https://doi.org/10.1101/217737",
    paper_year=2018,
    code_url="https://github.com/yanailab/knn-smoothing",
    image="openproblems-python-extras",
)
def knn_smoothing(adata, test=False):
    import knn_smooth

    adata.uns["method_code_version"] = check_version("knn_smooth")
    adata.obsm["train"] = (
        knn_smooth.knn_smoothing(np.array(adata.obsm["train"].transpose()), k=10)
    ).transpose()
    return adata

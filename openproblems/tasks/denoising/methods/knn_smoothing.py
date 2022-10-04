from ....tools.decorators import method
from ....tools.utils import check_version


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
    import numpy as np

    adata.uns["method_code_version"] = check_version("knn_smooth")
    X = adata.obsm["train"].transpose().toarray()
    X = X.astype(np.float64)
    adata.obsm["denoised"] = (knn_smooth.knn_smoothing(X, k=10)).transpose()
    return adata

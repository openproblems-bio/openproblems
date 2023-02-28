from ....tools.decorators import method
from ....tools.utils import check_version


@method(
    method_name="Iterative KNN smoothing",
    method_summary=(
        "Iterative kNN-smoothing is a method to repair or denoise noisy scRNA-seq"
        " expression matrices. Given a scRNA-seq expression matrix, KNN-smoothing first"
        " applies initial normalisation and smoothing. Then, a chosen number of"
        " principal components is used to calculate Euclidean distances between cells."
        " Minimally sized neighbourhoods are initially determined from these Euclidean"
        " distances, and expression profiles are shared between neighbouring cells."
        " Then, the resultant smoothed matrix is used as input to the next step of"
        " smoothing, where the size (k) of the considered neighbourhoods is increased,"
        " leading to greater smoothing. This process continues until a chosen maximum k"
        " value has been reached, at which point the iteratively smoothed object is"
        " then optionally scaled to yield a final result."
    ),
    paper_name=(
        "K-nearest neighbor smoothing for high-throughput single-cell RNA-Seq data"
    ),
    paper_reference="wagner2018knearest",
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

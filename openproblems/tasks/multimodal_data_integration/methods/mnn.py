import numpy as np
import scprep

from sklearn.decomposition import TruncatedSVD
from ....tools.normalize import log_cpm, log_scran_pooling
from ....tools.decorators import method
from ....tools.utils import check_version


def _mnn(adata, n_svd=100):
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri

    numpy2ri.activate()
    scprep.run.install_bioconductor("batchelor")

    if min(adata.X.shape) <= n_svd:
        n_svd = min(adata.X.shape) - 1
    if min(adata.obsm["mode2"].shape) <= n_svd:
        n_svd = min(adata.obsm["mode2"].shape) - 1
    X_pca = TruncatedSVD(n_svd).fit_transform(adata.X)
    Y_pca = TruncatedSVD(n_svd).fit_transform(adata.obsm["mode2"])

    ro.globalenv["expr"] = np.vstack([X_pca, Y_pca]).T
    ro.globalenv["batch"] = np.concatenate(
        [np.full(X_pca.shape[0], 1), np.full(Y_pca.shape[0], 2)]
    ).tolist()

    ro.r("batch <- as.integer(batch)")
    ro.r(
        "out <- as.matrix(SummarizedExperiment::assay(batchelor::fastMNN(expr, batch = batch), 'reconstructed'))"
    )

    XY_corrected = ro.globalenv["out"].T
    ro.r("rm(list=ls())")
    adata.obsm["aligned"] = XY_corrected[: X_pca.shape[0]]
    adata.obsm["mode2_aligned"] = XY_corrected[X_pca.shape[0] :]


@method(
    method_name="Mutual Nearest Neighbors (log CPM)",
    paper_name="Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors",
    paper_url="https://www.nature.com/articles/nbt.4091",
    paper_year=2018,
    code_url="https://github.com/LTLA/batchelor",
    code_version=check_version("rpy2"),
)
def mnn_log_cpm(adata, n_svd=100):
    log_cpm(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _mnn(adata, n_svd=n_svd)


@method(
    method_name="Mutual Nearest Neighbors (log scran)",
    paper_name="Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors",
    paper_url="https://www.nature.com/articles/nbt.4091",
    paper_year=2018,
    code_url="https://github.com/LTLA/batchelor",
    code_version=check_version("rpy2"),
)
def mnn_log_scran_pooling(adata, n_svd=100):
    log_scran_pooling(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    _mnn(adata, n_svd=n_svd)

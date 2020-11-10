import scprep

from ....tools.normalize import log_cpm, log_scran_pooling
from ....tools.decorators import method
from ....tools.utils import check_version


_mnn = scprep.run.RFunction(
    setup="""
    library(SingleCellExperiment)
    library(Matrix)
    library(sparsesvd)
    library(batchelor)
    """,
    args="sce, n_svd=100",
    body="""
    assay(sce, "X") <- as(assay(sce, "X"), "CsparseMatrix")
    reducedDim(sce, "mode2") <- as(reducedDim(sce, "mode2"), "CsparseMatrix")
    n_svd <- min(
      n_svd,
      dim(assay(sce, "X")),
      dim(reducedDim(sce, "mode2"))
    )
    X_svd <- sparsesvd(t(assay(sce, "X")), n_svd)
    X_svd <- X_svd[[2]] %*% diag(X_svd[[1]])
    Y_svd <- sparsesvd(reducedDim(sce, "mode2"), n_svd)
    Y_svd <- Y_svd[[2]] %*% diag(Y_svd[[1]])
    XY_svd <- rbind(X_svd, Y_svd)
    batch <- c(rep(1, nrow(X_svd)), rep(2, nrow(Y_svd)))
    XY_recons <- t(assay(fastMNN(t(XY_svd), batch = batch), 'reconstructed'))
    X_recons <- XY_recons[1:nrow(X_svd),]
    Y_recons <- XY_recons[(nrow(X_svd)+1):nrow(XY_svd),]
    reducedDim(sce, "aligned") <- as.matrix(X_recons)
    reducedDim(sce, "mode2_aligned") <- as.matrix(Y_recons)
    sce
    """,
)


@method(
    method_name="Mutual Nearest Neighbors (log CPM)",
    paper_name="Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors",
    paper_url="https://www.nature.com/articles/nbt.4091",
    paper_year=2018,
    code_url="https://github.com/LTLA/batchelor",
    code_version=check_version("rpy2"),
    image="openproblems-r-extras",
)
def mnn_log_cpm(adata, n_svd=100):
    log_cpm(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    return _mnn(adata, n_svd=n_svd)


@method(
    method_name="Mutual Nearest Neighbors (log scran)",
    paper_name="Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors",
    paper_url="https://www.nature.com/articles/nbt.4091",
    paper_year=2018,
    code_url="https://github.com/LTLA/batchelor",
    code_version=check_version("rpy2"),
    image="openproblems-r-extras",
)
def mnn_log_scran_pooling(adata, n_svd=100):
    log_scran_pooling(adata)
    log_cpm(adata, obsm="mode2", obs="mode2_obs", var="mode2_var")
    return _mnn(adata, n_svd=n_svd)

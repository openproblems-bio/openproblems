# Mutual Nearest Neighbors
# @param sce SingleCellExperiment
# @return sce SingleCellExperiment

# Libraries
library(SingleCellExperiment)
library(Matrix)
library(sparsesvd)
library(batchelor)

# Parameters
n_svd <- 100

# Run MNN
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
XY_recons <- t(assay(fastMNN(t(XY_svd), batch = batch), "reconstructed"))
X_recons <- XY_recons[1:nrow(X_svd), ]
Y_recons <- XY_recons[(nrow(X_svd) + 1):nrow(XY_svd), ]
reducedDim(sce, "aligned") <- as.matrix(X_recons)
reducedDim(sce, "mode2_aligned") <- as.matrix(Y_recons)

# Return
sce

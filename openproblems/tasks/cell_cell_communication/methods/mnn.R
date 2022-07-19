# Mutual Nearest Neighbors
# @param sce SingleCellExperiment
# @return sce SingleCellExperiment

# Libraries
library(SingleCellExperiment)
library(liana)

# Parameters
n_svd <- 100

# Convert data to friendly sparse format
assay(sce, "X") <- as(assay(sce, "X"), "CsparseMatrix")
reducedDim(sce, "mode2") <- as(reducedDim(sce, "mode2"), "CsparseMatrix")

# Check parameters
n_svd <- min(
  n_svd,
  dim(assay(sce, "X")),
  dim(reducedDim(sce, "mode2"))
)

# Run SVD
mode1_svd <- sparsesvd(t(assay(sce, "X")), n_svd)
mode1_svd <- mode1_svd[[2]] %*% diag(mode1_svd[[1]])
mode2_svd <- sparsesvd(reducedDim(sce, "mode2"), n_svd)
mode2_svd <- mode2_svd[[2]] %*% diag(mode2_svd[[1]])

# Run MNN
combined_svd <- rbind(mode1_svd, mode2_svd)
batch <- c(rep(1, nrow(mode1_svd)), rep(2, nrow(mode2_svd)))
sce_mnn <- fastMNN(t(combined_svd), batch = batch)

# Extract result
combined_recons <- t(assay(sce_mnn, "reconstructed"))
mode1_recons <- combined_recons[seq_len(nrow(mode1_svd)), ]
mode2_recons <- combined_recons[seq(nrow(mode1_svd) + 1, nrow(combined_svd)), ]

# Store in output SingleCellExperiment
reducedDim(sce, "aligned") <- as.matrix(mode1_recons)
reducedDim(sce, "mode2_aligned") <- as.matrix(mode2_recons)

# Return
sce

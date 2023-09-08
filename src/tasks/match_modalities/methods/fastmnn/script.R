library(anndata, warn.conflicts = FALSE)
library(Matrix, warn.conflicts = FALSE)
requireNamespace("batchelor", quietly = TRUE)

## VIASH START
par <- list(
  input_mod1 = "resources_test/common/multimodal/dataset_mod1.h5ad",
  input_mod2 = "resources_test/common/multimodal/dataset_mod2.h5ad",
  output_mod1 = "output_mod1.h5ad",
  output_mod2 = "output_mod2.h5ad"
)
## VIASH END

cat("Reading input h5ad file\n")
adata_mod1 <- read_h5ad(par$input_mod1)
adata_mod2 <- read_h5ad(par$input_mod2)

cat("Running MNN\n")
sce_mnn <- batchelor::fastMNN(
  t(adata_mod1$obsm[["X_svd"]]),
  t(adata_mod2$obsm[["X_svd"]])
)

cat("Storing output\n")
combined_recons <- t(SummarizedExperiment::assay(sce_mnn, "reconstructed"))
mode1_recons <- combined_recons[seq_len(nrow(adata_mod1$obsm[["X_svd"]])), , drop = FALSE]
mode2_recons <- combined_recons[-seq_len(nrow(adata_mod1$obsm[["X_svd"]])), , drop = FALSE]

adata_mod1$obsm[["integrated"]] <- as.matrix(mode1_recons)
adata_mod2$obsm[["integrated"]] <- as.matrix(mode2_recons)

cat("Writing to file\n")
adata_mod1$uns["method_id"] <- meta$functionality_name
adata_mod2$uns["method_id"] <- meta$functionality_name

yyy <- adata_mod1$write_h5ad(par$output_mod1, compression = "gzip")
zzz <- adata_mod2$write_h5ad(par$output_mod2, compression = "gzip")

## VIASH START
par <- list(
  input = "src/label_projection/resources/pancreas/toy_preprocessed_data.h5ad",
  output = "output.scran.h5ad"
)
## VIASH END

cat(">> Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
library(scran, warn.conflicts = FALSE)
library(BiocParallel, warn.conflicts = FALSE)

cat(">> Load data\n")
adata <- anndata::read_h5ad(par$input)
sce <- scran::calculateSumFactors(t(adata$X), min.mean=0.1, BPPARAM=BiocParallel::MulticoreParam())
adata$obs[["size_factors"]] <- sce
adata$X <- log1p(t(t(adata$X) * sce))

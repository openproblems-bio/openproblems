library(SingleCellExperiment)
library(batchelor)

expr <- assay(sce, "counts")

sce <- fastMNN(expr, batch = colData(sce)[[batch]], k = k, d = n_pca)

# return
if (return_features) {
  as(assay(sce, "reconstructed"), "dgCMatrix")
} else {
  reducedDim(sce, "corrected")
}

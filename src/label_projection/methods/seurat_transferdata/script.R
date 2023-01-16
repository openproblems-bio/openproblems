cat(">> Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
requireNamespace("Seurat", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)
library(magrittr, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input_train = "resources_test/label_projection/pancreas/train.h5ad",
  input_test = "resources_test/label_projection/pancreas/test.h5ad",
  output = "output.h5ad"
)
## VIASH END

cat(">> Load input data\n")
input_train <- read_h5ad(par$input_train)
input_test <- read_h5ad(par$input_test)

# sce_train <- zellkonverter::readH5AD(par$input_train)
# obj_train <- Seurat::as.Seurat(sce_train, data = "normalized")
# sce_test <- zellkonverter::readH5AD(par$input_test)
# obj_test <- Seurat::as.Seurat(sce_test, data = "normalized")

cat(">> Converting AnnData to Seurat\n")
anndataToSeurat <- function(adata) {
  # interpreted from https://github.com/satijalab/seurat/blob/v3.1.0/R/objects.R
  obj <-
    SeuratObject::CreateSeuratObject(
      counts = as(Matrix::t(adata$layers[["counts"]]), "CsparseMatrix")
    ) %>%
    SeuratObject::SetAssayData(
      slot = "data",
      new.data = as(Matrix::t(adata$layers[["normalized"]]), "CsparseMatrix")
    ) %>%
    SeuratObject::AddMetaData(
      adata$obs
    )

  # set hvg
  SeuratObject::VariableFeatures(obj) <- adata$var_names[adata$var[["hvg"]]]

  # set embedding
  # could add loadings and stdev
  embed <- SeuratObject::CreateDimReducObject(
    embeddings = adata$obsm[["X_pca"]],
    key = "PC_"
  )
  obj[["pca"]] <- embed

  # return
  obj
}

obj_train <- anndataToSeurat(input_train)
obj_test <- anndataToSeurat(input_test)

cat(">> Find transfer anchors\n")
npcs <- ncol(obj_train[["pca"]])
anchors <- Seurat::FindTransferAnchors(
  reference = obj_train,
  query = obj_test,
  npcs = npcs,
  dims = seq_len(npcs),
  verbose = FALSE
)

cat(">> Predict on test data\n")
query <- Seurat::TransferData(
  anchorset = anchors,
  reference = obj_train,
  query = obj_test,
  refdata = list(labels = "label"),
  verbose = FALSE
)
input_test$obs[["label_pred"]] <- query$predicted.labels[input_test$obs_names]

cat(">> Write output to file\n")
input_test$uns[["method_id"]] <- meta[["functionality_name"]]
input_test$write_h5ad(par$output, compression = "gzip")

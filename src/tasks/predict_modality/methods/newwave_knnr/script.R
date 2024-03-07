cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("NewWave", quietly = TRUE)
requireNamespace("FNN", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)

## VIASH START
path <- "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/"
par <- list(
  input_train_mod1 = paste0(path, "train_mod1.h5ad"),
  input_test_mod1 = paste0(path, "test_mod1.h5ad"),
  input_train_mod2 = paste0(path, "train_mod2.h5ad"),
  output = "output.h5ad",
  newwave_maxiter = 40L,
  newwave_ngene = 200L,
  newwave_ncell = 200L,
  n_neighbors = 20L
)
meta <- list(functionality_name = "foo")
## VIASH END

print(par)

n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)

method_id <- meta$functionality_name

cat("Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)

# fetch batch labels
batch1 <- c(as.character(input_train_mod1$obs$batch), as.character(input_test_mod1$obs$batch))
batch2 <- as.character(input_train_mod1$obs$batch)

# create SummarizedExperiment object
data1 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    counts = as(
      cbind(
        t(input_train_mod1$layers[["counts"]]),
        t(input_test_mod1$layers[["counts"]])
      ),
      "CsparseMatrix"
    )
  ),
  colData = data.frame(batch = factor(batch1))
)
data1 <- data1[Matrix::rowSums(SummarizedExperiment::assay(data1)) > 0, ]
rm(input_train_mod1, input_test_mod1)
gc()

cat("Running NewWave on mod1\n")
res1 <- NewWave::newWave(
  data1,
  X = "~batch",
  verbose = TRUE,
  K = 10,
  maxiter_optimize = par$newwave_maxiter,
  n_gene_par = min(par$newwave_ngene, nrow(data1)),
  n_cell_par = min(par$newwave_ncell, ncol(data1)),
  commondispersion = FALSE
)
dr_mod1 <- SingleCellExperiment::reducedDim(res1)
colnames(dr_mod1) <- paste0("comp_", seq_len(ncol(dr_mod1)))
rm(data1)
gc()

# split DR matrices
train_ix <- seq_along(batch2)
dr_mod1_train <- dr_mod1[train_ix, , drop = FALSE]
dr_mod1_test <- dr_mod1[-train_ix, , drop = FALSE]


cat("Predicting for each column in modality 2\n")
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)

# precompute knn indices
knn_ix <- FNN::get.knnx(
  dr_mod1_train,
  dr_mod1_test,
  k = min(nrow(dr_mod1_train), par$n_neighbors)
)$nn.index

# perform knn regression.
pred <- input_train_mod2$layers[["normalized"]][knn_ix[, 1], , drop = FALSE]
if (par$n_neighbors > 1) {
  for (k in seq(2, par$n_neighbors)) {
    pred <- pred + input_train_mod2$layers[["normalized"]][knn_ix[, k], , drop = FALSE]
  }
}
pred <- pred / par$n_neighbors
rownames(pred) <- rownames(dr_mod1_test)

cat("Creating outputs object\n")
out <- anndata::AnnData(
  layers = list(normalized = pred),
  shape = dim(pred),
  uns = list(
    dataset_id = input_train_mod2$uns[["dataset_id"]],
    method_id = meta$functionality_name
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

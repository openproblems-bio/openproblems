cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
path <- "output/datasets/predict_modality/openproblems_bmmc_multiome_phase1_mod1/openproblems_bmmc_multiome_phase1_mod1.censor_dataset.output_"
par <- list(
  input_train_mod1 = paste0(path, "train_mod1.h5ad"),
  input_test_mod1 = paste0(path, "test_mod1.h5ad"),
  input_train_mod2 = paste0(path, "train_mod2.h5ad"),
  output = "output.h5ad",
  n_pcs = 4L,
  n_neighbors = 3,
  distance_method = "pearson"
)
## VIASH END

cat("Reading mod1 h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
dataset_id <- input_train_mod1$uns[["dataset_id"]]

# subset to HVG to reduce memory consumption
train_mod1_sd <- proxyC::colSds(input_train_mod1$layers[["normalized"]])
ix <- order(train_mod1_sd, decreasing = TRUE)[seq_len(min(1000, length(train_mod1_sd)))]
input_train_mod1 <- input_train_mod1[,ix]$copy()
gc()

# subset to HVG to reduce memory consumption
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
input_test_mod1 <- input_test_mod1[,ix]$copy()
gc()

cat("Performing DR on the mod1 values\n")
# LMDS is more efficient than regular MDS because
# it does not compure a square distance matrix.
dr_mod1 <- lmds::lmds(
  rbind(input_train_mod1$layers[["normalized"]], input_test_mod1$layers[["normalized"]]),
  ndim = par$n_pcs,
  distance_method = par$distance_method
)

ix <- seq_len(nrow(input_train_mod1))
dr_mod1_train <- dr_mod1[ix, , drop = FALSE]
dr_mod1_test <- dr_mod1[-ix, , drop = FALSE]

# remove previous objects to save memory
rm(input_train_mod1, input_test_mod1)
gc()

cat("Reading mod2 h5ad files\n")
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)

cat("Predicting for each column in modality 2\n")
# precompute knn indices
knn_ix <- FNN::get.knnx(
  dr_mod1_train,
  dr_mod1_test,
  k = par$n_neighbors
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

out <- anndata::AnnData(
  layers = list(normalized = pred),
  shape = dim(pred),
  uns = list(
    dataset_id = dataset_id,
    method_id = meta$functionality_name
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
requireNamespace("pbapply", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
path <- "output/datasets/predict_modality/openproblems_bmmc_multiome_phase1_rna/openproblems_bmmc_multiome_phase1_rna.censor_dataset.output_"
par <- list(
  input_train_mod1 = paste0(path, "train_mod1.h5ad"),
  input_test_mod1 = paste0(path, "test_mod1.h5ad"),
  input_train_mod2 = paste0(path, "train_mod2.h5ad"),
  output = "output.h5ad",
  n_pcs = 20L,
  n_trees = 50L
)
meta <- list(functionality_name = "foo")
## VIASH END

n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)

cat("Reading mod1 files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)

dataset_id <- input_train_mod1$uns[["dataset_id"]]

cat("Performing DR on the mod1 values\n")
dr <- lmds::lmds(
  rbind(input_train_mod1$layers[["normalized"]], input_test_mod1$layers[["normalized"]]), 
  ndim = par$n_pcs,
  distance_method = par$distance_method
)

ix <- seq_len(nrow(input_train_mod1))
dr_train <- as.data.frame(dr[ix, , drop = FALSE])
dr_test <- as.data.frame(dr[-ix, , drop = FALSE])
dr_train <- dr[ix, , drop = FALSE]
dr_test <- dr[-ix, , drop = FALSE]

rm(input_train_mod1, input_test_mod1)
gc()


cat("Reading mod2 files\n")
X_mod2 <- anndata::read_h5ad(par$input_train_mod2)$layers[["normalized"]]

cat("Predicting for each column in modality 2\n")
preds <- pbapply::pblapply(
  seq_len(ncol(X_mod2)),
  cl = n_cores,
  function(i) {
    y <- X_mod2[, i]
    uy <- unique(y)
    if (length(uy) > 1) {
      rf <- ranger::ranger(
        x = dr_train,
        y = y,
        num.trees = par$n_trees
      )
      stats::predict(rf, dr_test)$prediction
    } else {
      rep(uy, nrow(dr_test))
    }
  }
)

cat("Creating outputs object\n")
prediction <- Matrix::Matrix(do.call(cbind, preds), sparse = TRUE)
rownames(prediction) <- rownames(dr_test)
colnames(prediction) <- colnames(X_mod2)

out <- anndata::AnnData(
  layers = list(normalized = prediction),
  shape = dim(prediction),
  uns = list(
    dataset_id = dataset_id,
    method_id = meta$functionality_name
  )
)


cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

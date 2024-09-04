cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_train_mod1 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/train_mod1.h5ad",
  input_test_mod1 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/test_mod1.h5ad",
  input_train_mod2 = "resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/train_mod2.h5ad",
  output = "output.h5ad"
)
meta <- list(functionality_name = "foo")
## VIASH END

cat("Reading h5ad files\n")
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)

cat("Creating outputs object\n")
sample_ix <- sample.int(nrow(input_train_mod2), nrow(input_test_mod1), replace = TRUE)
prediction <- input_train_mod2$layers[["normalized"]][sample_ix, , drop = FALSE]
rownames(prediction) <- rownames(input_test_mod1)

out <- anndata::AnnData(
  layers = list(normalized = prediction),
  shape = dim(prediction),
  uns = list(
    dataset_id = input_train_mod2$uns[["dataset_id"]],
    method_id = meta[["functionality_name"]]
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
